import argparse
import sys
import json
import shutil
import subprocess
import tempfile
import time
import os
from pathlib import Path
from .data_models import Candidate, read_json
from .sifi_adapter import import_sifi, provenance
from .efficiency_analysis import summarize
from .offtarget_reports import load_evaluations
from .reporting import write_reports
from .config import Settings, inventory, doctor
from .handoff import prepare, collect
from .automation import run_offtargets
from .reuse import offtarget_reuse_session
from .legacy_session import shared_sifi_session, shared_directory
from .window_explorer import explore_windows
from .window_reporting import write_window_outputs
from .automatic_window_selection import select_windows
from .automatic_window_reporting import write_selection_outputs
from .pipeline_workflow import run_integrated_workflow
from .console_output import print_workflow_summary
from .priority_scoring import priority_key, length_scores
from .legacy_runtime import detect_sifi_runtime, find_sifi_root


def protect_output(out: Path, database_dir: Path | None = None, sifi_root: Path | None = None):
    protected = []
    try:
        protected.append(find_sifi_root(sifi_root))
    except ValueError:
        pass
    if database_dir:
        protected.append(database_dir)
    for directory in protected:
        if out.resolve().is_relative_to(directory.resolve()):
            raise ValueError("La salida debe estar fuera del código original y de las bases")


def utility(command, argv):
    parser = argparse.ArgumentParser(prog=f"rnai-designer {command}")
    parser.add_argument("--config", type=Path)
    parser.add_argument("--db-dir", type=Path)
    parser.add_argument("--out", type=Path)
    if command == "doctor":
        parser.add_argument("--sifi-root", type=Path)
        parser.add_argument("--legacy-conda-env")
    if command == "databases":
        parser.add_argument("--inspect", action="store_true")
    if command == "prepare":
        parser.add_argument("--candidate", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        settings = Settings.load(args.config, args.db_dir)
        if command == "doctor":
            settings.sifi_root = args.sifi_root or settings.sifi_root
            settings.legacy_conda_env = args.legacy_conda_env or settings.legacy_conda_env
        if args.out:
            protect_output(args.out, settings.database_dir, settings.sifi_root)
        if command == "databases":
            result = inventory(settings, args.inspect)
        elif command == "doctor":
            result = doctor(settings)
        else:
            if args.out is None:
                raise ValueError("prepare requiere --out")
            result = prepare(Candidate.load(args.candidate), args.out)
        text = json.dumps(result, ensure_ascii=False, indent=2)
        if args.out and command != "prepare":
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(text, encoding="utf-8")
        print(text)
        return 0
    except (OSError, ValueError, TypeError, KeyError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2


def windows_command(argv):
    parser = argparse.ArgumentParser(
        prog="rnai-designer windows",
        description="Reconstruye regiones verdes de RNAi Designer y compara ventanas descriptivas")
    parser.add_argument("--report", type=Path, required=True, help="report.json integrado")
    parser.add_argument("--window-length", type=int, required=True, help="Longitud en nucleótidos")
    parser.add_argument("--step", type=int, default=1, help="Desplazamiento entre inicios de ventana")
    parser.add_argument("--sirna-size", type=int, help="Solo si no está en el manifiesto ni puede inferirse")
    parser.add_argument("--select-start", type=int, help="Inicio 1-based de la ventana a detallar")
    parser.add_argument("--out", type=Path, required=True, help="Directorio de salida nuevo o vacío")
    args = parser.parse_args(argv)
    try:
        protect_output(args.out)
        report = read_json(args.report)
        analysis = explore_windows(report, args.window_length, args.step,
                                   args.sirna_size, args.select_start)
        analysis = write_window_outputs(report, analysis, args.report, args.out)
    except (OSError, ValueError, TypeError, KeyError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    print(f'{len(analysis["green_regions"])} regiones verdes; '
          f'{len(analysis["windows"])} ventanas generadas en {args.out}.')
    if analysis["selected_window"]:
        selected = analysis["selected_window"]
        print(f'Detalle: {selected["window_id"]}; {selected["sirna_start_count"]} siRNAs; '
              f'{selected["offtarget_alignment_count"]} alineamientos off-target superpuestos.')
    for warning in analysis["warnings"][:10]:
        print(f"Advertencia: {warning}", file=sys.stderr)
    if len(analysis["warnings"]) > 10:
        print("Ver todas las advertencias en window_analysis.json", file=sys.stderr)
    return 0


def select_windows_command(argv):
    parser = argparse.ArgumentParser(
        prog="rnai-designer select-windows",
        description="Prioriza ventanas por etiquetas si-Fi, GC y diversidad; exige filtro off-target para liberar secuencia")
    parser.add_argument("--report", type=Path, required=True, help="report.json integrado")
    parser.add_argument("--min-length", type=int, default=180, help="Longitud mínima (default: 180 nt)")
    parser.add_argument("--max-length", type=int, default=250, help="Longitud máxima (default: 250 nt)")
    parser.add_argument("--sirna-size", type=int, help="Solo si no está en el manifiesto ni puede inferirse")
    parser.add_argument("--out", type=Path, required=True, help="Directorio de salida nuevo o vacío")
    parser.add_argument("--length-penalty-weight", type=float, default=30.0)
    args = parser.parse_args(argv)
    try:
        protect_output(args.out)
        report = read_json(args.report)
        result = select_windows(report, args.min_length, args.max_length, args.sirna_size, args.length_penalty_weight)
        winner = result["selected_window"]
        required_species = {"Homo sapiens", "Zea mays", "Apis mellifera"}
        by_name = {item.get("organism"): item for item in report.get("offtarget_evaluations", [])}
        passes = bool(winner) and set(by_name) == required_species and len(evaluations) == 3 and all(
            by_name[name].get("criterion", {}).get("status") == "criterion_met" and
            by_name[name].get("criterion", {}).get("exact_21nt_rows") == 0 for name in required_species)
        result["selected_window_offtarget_filter"] = {
            "passes": passes, "status": "criterion_met" if passes else "not_evaluated_or_failed",
            "scope": "El reporte del gen completo cubre los siRNAs de 21 nt de la ventana; requiere tres especies y cero hits exactos."}
        if winner:
            winner["sequence_release_status"] = "released_after_filter" if passes else "withheld_filter_failed_or_indeterminate"
            if not passes:
                winner.pop("sequence", None)
                winner["sequence_withheld"] = True
        result = write_selection_outputs(report, result, args.report, args.out)
    except (OSError, ValueError, TypeError, KeyError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    print(f'{result["evaluated_window_count"]} ventanas evaluadas en {len(result["green_regions"])} regiones verdes.')
    selected = result["selected_window"]
    if selected:
        print(f'Ventana priorizada: {selected["window_id"]}; prioridad heuristica {selected["efficacy_priority_score"]:.2f}/100; '
              f'filtro off-target aprobado: {result["selected_window_offtarget_filter"]["passes"]}.')
        if result["selected_window_offtarget_filter"]["passes"]:
            print(f'Secuencia aprobada por el criterio computacional: {selected["sequence"]}')
        else:
            print("Secuencia retenida: el reporte no acredita cero coincidencias exactas de 21 nt en las tres especies.")
    else:
        if result["selection_status"] == "no_evaluable_efficacy_metrics":
            print("Hay ventanas, pero faltan clasificaciones si-Fi evaluables para priorizarlas.")
        else:
            print("No hay ventanas dentro del rango y las regiones verdes disponibles.")
    for warning in result["warnings"][:10]:
        print(f"Advertencia: {warning}", file=sys.stderr)
    return 0


def _legacy_sifi_context(args, candidate, out, settings):
    sifi_root = find_sifi_root(args.sifi_root or getattr(settings, "sifi_root", None))
    args.sifi_root = sifi_root
    protect_output(out, sifi_root)
    runtime = getattr(args, "_sifi_runtime", None)
    if runtime is None:
        runtime = detect_sifi_runtime(sifi_root, args.legacy_conda_env or getattr(settings, "legacy_conda_env", None))
        args._sifi_runtime = runtime
    if not runtime.get("available"):
        raise ValueError("No se pudo detectar un entorno Conda ejecutable para si-Fi: " + runtime.get("reason", "sin detalles"))
    for component, path in (("Bowtie", settings.bowtie), ("RNAplfold", settings.rnaplfold)):
        if not path or not path.is_file():
            raise ValueError(component + " no encontrado; instalarlo o configurar su ejecutable.")
    if not settings.database_dir.is_dir():
        raise ValueError("Bases de datos no encontradas: " + str(settings.database_dir))
    query = out / "designer_handoff" / "query.fasta"
    if not query.is_file():
        raise ValueError("Falta la query preparada por workflow: " + str(query))
    inputs = out / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    source_path, mask_path = inputs / "designer_original.json", inputs / "green_mask.json"
    bridge = Path(__file__).resolve().parent / "legacy_bridge" / "run_sifi_headless.py"
    common = [str(runtime["conda"]), "run", "-p", runtime["prefix"], "python", str(bridge),
              "--sifi-root", str(sifi_root), "--query-fasta", str(query),
              "--database-prefix", settings.databases["Dalbulus maidis"]["prefix"],
              "--database-dir", str(settings.database_dir), "--bowtie", str(settings.bowtie),
              "--rnaplfold-dir", str(settings.rnaplfold.parent), "--out-json", str(source_path),
              "--out-mask", str(mask_path), "--sirna-size", str(args.legacy_sirna_size),
              "--mismatches", str(args.sifi_mismatches), "--accessibility-window", str(args.accessibility_window),
              "--end-stability-threshold", str(args.end_stability_threshold),
              "--target-accessibility-threshold", str(args.target_accessibility_threshold)]
    if shared_directory.get() is not None:
        common += ["--shared-dir", shared_directory.get()]
    return runtime, query, source_path, mask_path, common


def _list_legacy_targets(args, candidate, out, settings):
    runtime, query, source_path, mask_path, common = _legacy_sifi_context(args, candidate, out, settings)
    process = subprocess.run(common + ["--list-targets"], capture_output=True, text=True,
                             errors="replace", timeout=args.sifi_timeout,
                             creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    if process.returncode:
        raise ValueError("Bowtie no pudo enumerar main targets: " + (process.stderr or process.stdout)[-2500:])
    payload = next((json.loads(line) for line in reversed(process.stdout.splitlines())
                    if line.lstrip().startswith("{")), None)
    if not isinstance(payload, dict):
        raise ValueError("El adaptador si-Fi no devolvio la lista JSON de targets.")
    return payload


def _run_legacy_sifi(args, candidate, out, settings, main_target):
    started = time.perf_counter()
    runtime, query, source_path, mask_path, common = _legacy_sifi_context(args, candidate, out, settings)
    alignment_path = Path(str(source_path) + '.bowtie.tmp')
    cache_artifacts = {'designer_original.json': source_path, 'green_mask.json': mask_path,
                       'designer_original.json.bowtie.tmp': alignment_path}
    cache = None
    cache_hit = False
    loaded_records = loaded_mask = None
    if getattr(args, 'cache_dir', None):
        if args.cache_dir.resolve().is_relative_to(settings.database_dir.resolve()):
            raise ValueError('La caché debe estar fuera del directorio de bases originales')
        from .cache import ArtifactCache, content_key, file_digest
        from .config import index_files
        prefix = Path(runtime['prefix'])
        interpreter = prefix / ('python.exe' if os.name == 'nt' else 'bin/python')
        files = [interpreter, Path(runtime['conda']),
                 Path(__file__).parent/'legacy_bridge'/'run_sifi_headless.py',
                 *args.sifi_root.rglob('*.py'), *prefix.glob('conda-meta/*.json'),
                 *index_files(settings.database_dir, settings.databases['Dalbulus maidis']['prefix'])]
        for tool in (settings.bowtie, settings.rnaplfold):
            files.extend(p for p in tool.parent.iterdir() if p.is_file() and p.suffix.lower() in ('.exe', '.dll'))
            files.append(tool)
        from . import __version__
        parameters = {'candidate': candidate.original, 'runtime': runtime, 'query_sha256': file_digest(query),
                      'rnai_designer_version': __version__,
                      'database_binding': settings.databases['Dalbulus maidis'],
                      'main_target': main_target, 'sirna_size': args.legacy_sirna_size,
                      'mismatches': args.sifi_mismatches, 'accessibility_window': args.accessibility_window,
                      'end_stability_threshold': args.end_stability_threshold,
                      'target_accessibility_threshold': args.target_accessibility_threshold}
        cache = ArtifactCache(args.cache_dir, content_key(parameters, set(files)))
        cache_hit = cache.restore(cache_artifacts)
    print(f'[si-Fi] {main_target}: ' + ('caché verificada' if cache_hit else 'ejecutando Conda'), file=sys.stderr, flush=True)
    if not cache_hit:
        from .atomic_io import temporary_directory
        with temporary_directory(source_path.parent, '.sifi.') as temporary:
            temporary_source = Path(temporary)/source_path.name
            temporary_mask = Path(temporary)/mask_path.name
            command = list(common)
            command[command.index('--out-json')+1] = str(temporary_source)
            command[command.index('--out-mask')+1] = str(temporary_mask)
            process = subprocess.run(command + ["--main-target", main_target], capture_output=True, text=True,
                                     errors="replace", timeout=args.sifi_timeout,
                                     creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            if process.returncode:
                raise ValueError("Ejecucion headless para {0} fallo: {1}".format(
                    main_target, (process.stderr or process.stdout)[-3000:]))
            loaded_records, loaded_mask = read_json(temporary_source), read_json(temporary_mask)
            if not isinstance(loaded_records, list) or not isinstance(loaded_mask, dict):
                raise ValueError('si-Fi no generó artefactos JSON completos')
            temporary_alignments = Path(str(temporary_source)+'.bowtie.tmp')
            if not temporary_alignments.is_file():
                raise ValueError('si-Fi no conservó la evidencia de alineamientos Bowtie')
            os.replace(temporary_source, source_path)
            os.replace(temporary_mask, mask_path)
            os.replace(temporary_alignments, alignment_path)
    query_id = query.read_text(encoding="utf-8").splitlines()[0].lstrip("> ").split()[0]
    metadata = {"schema_version": "1.0", "engine": "si-Fi original Python 2 headless adapter",
                "query_id": query_id, "gene_id": candidate.gene_id, "main_target": main_target,
                "runtime": {key: runtime.get(key) for key in ("name", "prefix", "python", "bits")},
                "parameters": {"sirna_size": args.legacy_sirna_size, "mismatches": args.sifi_mismatches,
                               "accessibility_window": args.accessibility_window,
                               "end_stability_threshold": args.end_stability_threshold,
                               "target_accessibility_threshold": args.target_accessibility_threshold,
                               "database_prefix": settings.databases["Dalbulus maidis"]["prefix"],
                               "main_targets": [main_target]},
                "tools": {"bowtie": str(settings.bowtie), "rnaplfold": str(settings.rnaplfold)},
                "query_fasta": str(query), "result_json": "inputs/designer_original.json",
                "green_mask": "inputs/green_mask.json", "headless": True,
                "bowtie_alignments": "inputs/designer_original.json.bowtie.tmp"}
    run_path = out / "designer_handoff" / "run.json"
    if run_path.exists():
        existing = read_json(run_path)
        if existing.get("status") not in ("not_evaluated", None):
            raise ValueError("Ya existe un run.json completado en " + str(out))
    from .atomic_io import write_json
    write_json(run_path, metadata)
    imported = import_sifi(source_path, candidate, metadata, _rows=loaded_records)
    imported["green_mask"] = read_json(mask_path) if loaded_mask is None else loaded_mask
    imported["green_mask_provenance"] = provenance(mask_path)
    imported["warnings"].extend(imported["green_mask"].get("warnings", []))
    if cache is not None and not cache_hit:
        cache.store(cache_artifacts)
    write_json(out/'headless_performance.json', {'elapsed_s': time.perf_counter()-started, 'cache_hit': cache_hit})
    return source_path, imported, metadata


def _safe_reference_slug(name):
    import re
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._") or "target"


@offtarget_reuse_session()
def _run_multi_reference_workflow(args, candidate, settings, out, main_targets):
    import csv
    import tempfile
    import shutil
    from .atomic_io import atomic_text, write_json
    from copy import copy
    references_root = out / "references"
    references_root.mkdir(parents=True, exist_ok=True)
    reference_rows = []
    combined_count = 0
    combined_fields = None
    combined_spool = tempfile.TemporaryFile(mode='w+', encoding='utf-8', newline='')
    for index, target in enumerate(dict.fromkeys(main_targets), 1):
        slug = "{0:02d}_{1}".format(index, _safe_reference_slug(target))
        reference_out = references_root / slug
        entry = {"main_target": target, "reference_id": target, "output_directory": str(reference_out.relative_to(out))}
        spool_position = combined_spool.tell()
        previous_fields, previous_count = combined_fields, combined_count
        try:
            reference_out.mkdir(parents=True, exist_ok=False)
            (reference_out / "candidate.json").write_text(
                json.dumps(candidate.original, ensure_ascii=False, indent=2), encoding="utf-8")
            prepare(candidate, reference_out / "designer_handoff")
            target_args = copy(args)
            target_args._sifi_runtime = args._sifi_runtime
            source_path, imported, metadata = _run_legacy_sifi(target_args, candidate, reference_out,
                                                                settings, target)
            run_path = reference_out / "designer_handoff" / "run.json"
            imported["candidate_source"] = provenance(reference_out / "candidate.json")
            result = run_integrated_workflow(candidate, imported, settings, reference_out, source_path,
                                             run_path, args.timeout, args.min_window, args.max_window,
                                             args.sirna_size, **({"length_penalty_weight": args.length_penalty_weight} if hasattr(args, "length_penalty_weight") else {}))
            summary = result["summary"]
            entry.update(status="complete", green_regions=summary["green_regions"],
                         evaluated_window_count=summary["evaluated_window_count"],
                         selected_window=summary["selected_window"],
                         selected_window_offtarget_filter=summary["selected_window_offtarget_filter"],
                         selection_method=result["selection"]["selection_method"],
                         warnings=summary["warnings"],
                         selection_status=summary.get("selection_status"),
                         offtarget_selection=summary.get("offtarget_selection"),
                         selected_window_offtarget_evaluations=summary.get("selected_window_offtarget_evaluations", []))
            windows = result['selection']['windows']
            if windows:
                if hasattr(windows, 'scalar_rows'):
                    fields = ['main_target', 'window_id', *windows.columns]
                    rows = ([target, *values] for values in windows.scalar_rows())
                else:
                    fields = ['main_target', *windows[0]]
                    rows = ([target, *(row.get(key) for key in fields[1:])] for row in windows)
                if combined_fields is None:
                    combined_fields = fields
                    csv.writer(combined_spool).writerow(fields)
                if fields != combined_fields:
                    raise ValueError('Esquemas de ventanas incompatibles entre referencias')
                csv.writer(combined_spool).writerows(rows)
                combined_count += len(windows)
            # Keep only the per-reference summary, not all windows of all references.
            del result
        except (OSError, ValueError, TypeError, KeyError, subprocess.TimeoutExpired) as error:
            # Discard a partially written reference; its prefix must never appear
            # as valid rows in the consolidated CSV after a recoverable failure.
            combined_spool.seek(spool_position)
            combined_spool.truncate()
            combined_fields, combined_count = previous_fields, previous_count
            entry.update(status="error", error=str(error), green_regions=[],
                         evaluated_window_count=0, selected_window=None, warnings=[str(error)])
            reference_out.mkdir(parents=True, exist_ok=True)
            (reference_out / "reference_error.json").write_text(
                json.dumps(entry, ensure_ascii=False, indent=2), encoding="utf-8")
        reference_rows.append(entry)

    completed_count = sum(row["status"] == "complete" for row in reference_rows)
    failed_count = sum(row["status"] == "error" for row in reference_rows)
    report = {"schema_version": "1.0", "candidate_id": candidate.gene_id,
              "status": ("complete" if completed_count and not failed_count else
                         "partial" if completed_count else "no_reference_completed"),
              "reference_count": len(reference_rows), "completed_reference_count": completed_count,
              "failed_reference_count": failed_count,
              "references": reference_rows,
              "selection_method": next((row["selection_method"] for row in reference_rows
                                        if row["status"] == "complete"),
                                       {"score_name": "efficacy_priority_score",
                                        "interpretation": "Heuristica exploratoria; no equivale a eficacia biologica validada."}),
              "combined_window_count": combined_count,
              "combined_windows_csv": "all_references_windows.csv", "table_schema_version": "2.0",
              "coordinate_note": "Las coordenadas son relativas al query de cada referencia y solo se comparan dentro del mismo main_target."}
    eligible = [r for r in reference_rows if r["status"] == "complete" and r.get("selected_window")
                and r["selected_window"].get("adjusted_priority_score", r["selected_window"].get("efficacy_priority_score")) is not None]
    best = min(eligible, key=lambda r: (priority_key(r["selected_window"]), r["main_target"])) if eligible else None
    report["selected_window"] = ({**best["selected_window"], "main_target": best["main_target"]} if best else None)
    report["selected_reference"] = best["main_target"] if best else None
    csv_path = out / "all_references_windows.csv"
    with atomic_text(csv_path, newline='') as handle:
        if combined_fields is None:
            csv.writer(handle).writerow(['main_target', 'region_id', 'start_1based', 'end_1based_inclusive', 'length_nt', 'window_class'])
        else:
            combined_spool.seek(0)
            shutil.copyfileobj(combined_spool, handle)
    combined_spool.close()
    write_json(out / 'reference_comparison.json', report)
    write_json(out / 'pipeline_result.json', report)
    return report

@shared_sifi_session()
def workflow_command(argv, *, input_metadata=None):
    parser = argparse.ArgumentParser(
        prog="rnai-designer workflow",
        description="Prepara un gen; ejecuta cada main target de si-Fi por separado y consolida ventanas/off-targets")
    parser.add_argument("--resume", type=Path, help="Carpeta creada por una preparación workflow previa")
    parser.add_argument("--gene-id", help="Identificador del gen candidato")
    sequence = parser.add_mutually_exclusive_group()
    sequence.add_argument("--sequence", help="Secuencia IUPAC o FASTA de una sola entrada")
    sequence.add_argument("--sequence-file", type=Path, help="Archivo FASTA de una sola entrada")
    parser.add_argument("--reference-source", help="Procedencia de la secuencia, si está disponible")
    parser.add_argument("--reference-version", help="Versión/accesión, si está disponible")
    parser.add_argument("--homolog-id", help="ID de homólogo upstream, si está disponible")
    parser.add_argument("--priority-score", type=float, help="Puntaje upstream ya calculado, si existe")
    parser.add_argument("--out", type=Path, help="Directorio nuevo para iniciar; en resume se usa --resume")
    parser.add_argument("--sifi", help="JSON de si-Fi o 'auto' para recuperar el temporal compatible")
    parser.add_argument("--run", type=Path, help="Manifiesto si-Fi completado; por omisión usa el de designer_handoff")
    parser.add_argument("--green-mask", type=Path, help="Sidecar exacto exportado desde general_helpers.get_target_data")
    parser.add_argument("--legacy-python", type=Path, help="Intérprete Python 2.7 con dependencias de si-Fi para exportar la máscara directamente")
    parser.add_argument("--run-sifi-headless", action="store_true", help="Ejecuta el cálculo original de si-Fi en Conda sin abrir la GUI")
    parser.add_argument("--legacy-conda-env", help="Entorno Conda; por omisión se detecta uno con 'sifi' en el nombre")
    parser.add_argument("--main-target", action="append", default=[], help="ID Bowtie de target principal; repetible")
    parser.add_argument("--list-targets-only", action="store_true", help="Enumera targets Bowtie y termina")
    parser.add_argument("--legacy-sirna-size", type=int, default=21)
    parser.add_argument("--sifi-mismatches", type=int, default=0)
    parser.add_argument("--accessibility-window", type=int, default=8)
    parser.add_argument("--end-stability-threshold", type=float, default=1.0)
    parser.add_argument("--target-accessibility-threshold", type=float, default=0.1)
    parser.add_argument("--sifi-timeout", type=float, default=1800)
    parser.add_argument("--sifi-root", type=Path, default=None,
                        help="Carpeta original si-Fi usada por --legacy-python")
    parser.add_argument("--search-dir", type=Path, action="append", help="Directorio temporal adicional para --sifi auto")
    parser.add_argument("--wait-for-sifi", action="store_true", help="Espera mientras se ejecuta si-Fi GUI y luego continúa automáticamente")
    parser.add_argument("--wait-timeout", type=float, default=3600, help="Límite de espera en segundos; default 3600")
    parser.add_argument("--config", type=Path)
    parser.add_argument("--db-dir", type=Path)
    parser.add_argument("--timeout", type=float, default=300, help="Timeout Bowtie por organismo, segundos")
    parser.add_argument("--sirna-size", type=int, help="Tamaño siRNA si el manifiesto no lo registra")
    parser.add_argument("--length-penalty-weight", type=float, default=30.0)
    parser.add_argument("--min-window", type=int, default=180)
    parser.add_argument("--max-window", type=int, default=250)
    parser.add_argument('--cache-dir', type=Path, help='Caché opcional de si-Fi headless verificada por contenido, parámetros, herramientas y entorno')
    parser.add_argument("--verbose-report", action="store_true", help="Mostrar el reporte detallado de consola")
    args = parser.parse_args(argv)
    try:
        length_scores(None, 180, args.length_penalty_weight)
    except ValueError as error:
        parser.error(str(error))

    try:
        settings = Settings.load(args.config, args.db_dir)
        if args.wait_timeout <= 0 or args.timeout <= 0 or args.sifi_timeout <= 0:
            raise ValueError("Los timeout deben ser positivos")
        if args.run_sifi_headless and (args.sifi or args.wait_for_sifi):
            raise ValueError("--run-sifi-headless no se combina con --sifi ni --wait-for-sifi")
        if args.run_sifi_headless and (args.green_mask or args.legacy_python):
            raise ValueError("La ejecución headless ya exporta su máscara; no use --green-mask/--legacy-python")
        if args.legacy_sirna_size < 1 or args.sifi_mismatches < 0:
            raise ValueError("--legacy-sirna-size debe ser positivo y --sifi-mismatches no negativo")
        if args.resume:
            out = args.resume.resolve()
            protect_output(out, settings.database_dir, args.sifi_root or getattr(settings, "sifi_root", None))
            candidate_path = out / "candidate.json"
            candidate = Candidate.load(candidate_path)
            run_path = args.run or (out / "designer_handoff" / "run.json")
            if not run_path.is_file():
                run_path = None
        else:
            if args.out is None or not args.gene_id or (args.sequence is None and args.sequence_file is None):
                raise ValueError("Para iniciar se requieren --gene-id, --sequence o --sequence-file y --out; o usar --resume")
            out = args.out.resolve()
            protect_output(out, settings.database_dir, args.sifi_root or getattr(settings, "sifi_root", None))
            if out.exists() and any(out.iterdir()):
                raise ValueError("El directorio de workflow ya tiene archivos; usar --resume o una carpeta nueva")
            raw_sequence = args.sequence if args.sequence is not None else args.sequence_file.read_text(encoding="utf-8-sig")
            candidate = Candidate.from_fields(args.gene_id, raw_sequence, args.reference_source,
                                               args.reference_version, args.homolog_id, args.priority_score)
            if input_metadata is not None:
                candidate.original["metadata"].update(input_metadata)
            out.mkdir(parents=True, exist_ok=True)
            candidate_path = out / "candidate.json"
            candidate_path.write_text(json.dumps(candidate.original, ensure_ascii=False, indent=2), encoding="utf-8")
            prepare(candidate, out / "designer_handoff")
            run_path = args.run or (out / "designer_handoff" / "run.json")

        if run_path is not None and not run_path.is_file():
            raise ValueError(f"No se encuentra el manifiesto de ejecución: {run_path}")
        run_metadata = read_json(run_path) if run_path else None
        if run_metadata is not None and not isinstance(run_metadata, dict):
            raise ValueError("El manifiesto --run debe ser objeto")
        if args.green_mask and args.legacy_python:
            raise ValueError("Elegí --green-mask o --legacy-python, no ambos")
        if not args.sifi and not args.wait_for_sifi and not args.run_sifi_headless:
            print(f'Preparación lista en {out}. Abrí {out / "designer_handoff" / "query.fasta"} en RNAi Designer, '
                  "ejecutá RNAi design contra la base Dalbulus maidis y completá designer_handoff/run.json.")
            print(f'Después reanudá con: python -m rnai_designer workflow --resume "{out}" --sifi auto')
            return 0

        settings = Settings.load(args.config, args.db_dir)
        if args.run_sifi_headless:
            protect_output(out, settings.database_dir, args.sifi_root or getattr(settings, "sifi_root", None))
            _legacy_sifi_context(args, candidate, out, settings)  # Detect once and reuse this runtime for all targets.
            target_payload = None
            if args.list_targets_only or not args.main_target:
                target_payload = _list_legacy_targets(args, candidate, out, settings)
                detected_targets = [entry["hit_name"] for entry in target_payload.get("available_targets", [])]
            else:
                detected_targets = []
            main_targets = (detected_targets if args.list_targets_only or not args.main_target
                            else list(dict.fromkeys(args.main_target)))
            if not args.list_targets_only and not args.main_target:
                print("Targets si-Fi detectados: " + (", ".join(main_targets) if main_targets else "ninguno"))
            if args.list_targets_only:
                print(json.dumps(target_payload, ensure_ascii=False, indent=2))
                return 0
            if not main_targets:
                raise ValueError("Bowtie no encontro main targets para esta secuencia.")
            for subdir in (out / "references", out / "offtarget_runs", out / "offtarget_selected_window",
                           out / "integrated", out / "window_selection"):
                if subdir.exists() and any(subdir.iterdir()):
                    raise ValueError("La carpeta de salida ya contiene resultados; usar una carpeta nueva")
            for artifact in (out / "reference_comparison.json", out / "all_references_windows.csv",
                             out / "pipeline_result.json"):
                if artifact.exists():
                    raise ValueError("Ya existe un resultado consolidado: " + str(artifact))
            consolidated = _run_multi_reference_workflow(args, candidate, settings, out, main_targets)
            print_workflow_summary(consolidated, out, candidate.sequence, verbose=args.verbose_report)
            return 0 if (consolidated["completed_reference_count"] > 0 and
                         consolidated["failed_reference_count"] == 0) else 2
        source_path = None
        search_dirs = args.search_dir or [Path(tempfile.gettempdir())]
        if args.sifi and args.sifi.lower() != "auto":
            source_path = Path(args.sifi).resolve()
            imported = import_sifi(source_path, candidate, run_metadata)
        else:
            deadline = time.monotonic() + (args.wait_timeout if args.wait_for_sifi else 0)
            if args.wait_for_sifi:
                print(f'Esperando el JSON temporal de si-Fi compatible con {candidate.gene_id}. '
                      "Ejecutá la GUI si-Fi y completá el análisis; Ctrl+C cancela la espera.")
            while True:
                try:
                    source_path, imported = collect(candidate, run_metadata, search_dirs)
                    break
                except ValueError as error:
                    if not args.wait_for_sifi or time.monotonic() >= deadline:
                        raise ValueError(f"No se encontraron resultados si-Fi compatibles: {error}")
                    time.sleep(min(3, max(.1, deadline - time.monotonic())))
            if args.wait_for_sifi and run_path is not None:
                run_metadata = read_json(run_path)
                if not isinstance(run_metadata, dict):
                    raise ValueError("El manifiesto --run debe ser objeto")
                imported = import_sifi(source_path, candidate, run_metadata)

        if args.green_mask:
            mask_path = args.green_mask.resolve()
            imported["green_mask"] = read_json(mask_path)
            imported["green_mask_provenance"] = provenance(mask_path)
            inputs = out / "inputs"
            inputs.mkdir(parents=True, exist_ok=True)
            mask_copy = inputs / "green_mask.json"
            if mask_path != mask_copy.resolve():
                source_bytes = mask_path.read_bytes()
                if mask_copy.exists() and mask_copy.read_bytes() != source_bytes:
                    raise ValueError("Ya existe una máscara distinta en inputs; usar una carpeta nueva")
                if not mask_copy.exists():
                    mask_copy.write_bytes(source_bytes)
        elif args.legacy_python:
            run_metadata = run_metadata or {}
            parameters = run_metadata.get("parameters") if isinstance(run_metadata.get("parameters"), dict) else {}
            mask_size = args.sirna_size or parameters.get("sirna_size")
            if mask_size is None:
                sizes = {len(row["sirna_sequence"]) for row in imported["records"]
                         if isinstance(row.get("sirna_sequence"), str)}
                if len(sizes) == 1:
                    mask_size = next(iter(sizes))
            if type(mask_size) is not int or mask_size < 1:
                raise ValueError("Para --legacy-python indicar --sirna-size o registrar sirna_size en run.json")
            if source_path is None:
                raise ValueError("No se conoce el archivo JSON de si-Fi para exportar la máscara original")
            mask_path = out / "inputs" / "green_mask.json"
            mask_path.parent.mkdir(parents=True, exist_ok=True)
            bridge = Path(__file__).resolve().parent / "legacy_bridge" / "export_green_mask.py"
            process = subprocess.run(
                [str(args.legacy_python), str(bridge), "--sifi-json", str(source_path),
                 "--sirna-size", str(mask_size), "--sifi-root", str(find_sifi_root(args.sifi_root or getattr(settings, "sifi_root", None))),
                 "--out", str(mask_path)], capture_output=True, text=True, errors="replace",
                timeout=120, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            if process.returncode != 0:
                raise ValueError("El puente Python 2 de si-Fi falló: " + (process.stderr or process.stdout)[-2000:])
            imported["green_mask"] = read_json(mask_path)
            imported["green_mask_provenance"] = provenance(mask_path)

        protect_output(out, settings.database_dir, args.sifi_root or getattr(settings, "sifi_root", None))
        for subdir in (out / "offtarget_runs", out / "offtarget_selected_window",
                       out / "integrated", out / "window_selection"):
            if subdir.exists() and any(subdir.iterdir()):
                raise ValueError(f"La carpeta de salida ya contiene archivos: {subdir}; iniciar otro workflow")
        imported["candidate_source"] = provenance(candidate_path)
        result = run_integrated_workflow(candidate, imported, settings, out, source_path, run_path,
                                         args.timeout, args.min_window, args.max_window, args.sirna_size,
                                         length_penalty_weight=args.length_penalty_weight)
    except (OSError, ValueError, TypeError, KeyError, subprocess.TimeoutExpired) as error:
        print(f"Error: {error}", file=sys.stderr)
        print("Workflow fallido o incompleto: sin ganador final. Motivo: " + str(error), file=sys.stderr)
        return 2

    summary = result["summary"]
    print_workflow_summary(summary, out, candidate.sequence, verbose=args.verbose_report)
    safety = summary["offtarget_21nt_filter"]
    selected_errors = any(entry["status"] == "error"
                          for entry in summary["selected_window_offtarget_evaluations"])
    return 2 if selected_errors or any(entry["run_status"] == "error" for entry in safety["by_species"]) else 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "batch":
        from .batch import batch_command
        return batch_command(argv[1:], workflow_command)
    if argv and argv[0] == "windows":
        return windows_command(argv[1:])
    if argv and argv[0] == "select-windows":
        return select_windows_command(argv[1:])
    if argv and argv[0] == "workflow":
        return workflow_command(argv[1:])
    if argv and argv[0] in ("databases", "doctor", "prepare"):
        return utility(argv[0], argv[1:])
    if argv and not argv[0].startswith("-") and argv[0] not in ("analyze", "collect"):
        from .simple_run import simple_command
        return simple_command(argv, workflow_command)
    if argv and argv[0] == "analyze":
        argv = ["--automate-offtargets", *argv[1:]]
    elif argv and argv[0] == "collect":
        argv = ["--sifi", "auto", *argv[1:]]
    parser = argparse.ArgumentParser(description="Importa si-Fi y automatiza búsquedas con su motor Bowtie")
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--sifi", type=Path, help="JSON original o auto para buscar por ID y secuencia")
    parser.add_argument("--run", type=Path, help="Manifiesto de procedencia de la ejecución primaria")
    parser.add_argument("--offtargets", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--db-dir", type=Path)
    parser.add_argument("--search-dir", type=Path, action="append")
    parser.add_argument("--automate-offtargets", action="store_true")
    parser.add_argument("--timeout", type=float, default=300)
    args = parser.parse_args(argv)
    if args.sifi is None and not args.automate_offtargets:
        parser.error("Se requiere --sifi o el comando analyze")
    if args.offtargets and args.automate_offtargets:
        parser.error("Elegí --offtargets para importar o analyze para ejecutar")
    if args.timeout <= 0:
        parser.error("--timeout debe ser positivo")
    try:
        settings = Settings.load(args.config, args.db_dir) if args.automate_offtargets or args.config or args.db_dir else None
        protect_output(args.out, settings.database_dir if settings else None)
        candidate = Candidate.load(args.candidate)
        run = read_json(args.run) if args.run else None
        if run is not None and not isinstance(run, dict):
            raise ValueError("El manifiesto --run debe ser objeto")
        if args.sifi == Path("auto"):
            source, imported = collect(candidate, run, args.search_dir)
            destination = args.out / "inputs" / "designer_original.json"
            destination.parent.mkdir(parents=True, exist_ok=True)
            if source.resolve() != destination.resolve():
                shutil.copyfile(source, destination)
            imported["recovered_copy"] = "inputs/designer_original.json"
        elif args.sifi:
            imported = import_sifi(args.sifi, candidate, run)
        else:
            imported = {"source": {}, "run": None, "records": [], "status": "not_evaluated",
                        "warnings": ["RNAi Designer no ejecutado; eficiencia y accesibilidad no evaluadas."]}
            if candidate.original["metadata"].get("simulated"):
                imported["warnings"].append("DATOS FICTICIOS: query artificial de demostración.")
            prepare(candidate, args.out / "designer_handoff")
        imported["candidate_source"] = provenance(args.candidate)
        if args.run:
            imported["run_source"] = provenance(args.run)
        summary = summarize(imported["records"])
        evaluations = (run_offtargets(candidate, settings, args.out / "offtarget_runs", args.timeout)
                       if args.automate_offtargets else load_evaluations(args.offtargets, candidate, tolerant=True))
        report = write_reports(args.out, candidate, imported, summary, evaluations)
    except (OSError, ValueError, TypeError, KeyError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    print(f'Reportes en {args.out}; {len(report["warnings"])} advertencias. '
          f'Criterio 21nt: {report["offtarget_criterion"]["label"]}. Seguridad: indeterminada.')
    for warning in report["warnings"][:10]:
        print(f"Advertencia: {warning}", file=sys.stderr)
    if len(report["warnings"]) > 10:
        print("Ver todas las advertencias en report.json", file=sys.stderr)
    return 2 if any(e["status"] == "error" for e in evaluations) else 0


if __name__ == "__main__":
    sys.exit(main())
