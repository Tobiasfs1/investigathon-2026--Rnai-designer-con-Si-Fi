"""Render documentation from recorded measurements, never from estimated results."""
import json
import hashlib
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/'benchmarks'/'results'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def replace_section(text, name, content):
    start, end = '<!-- '+name+'_START -->', '<!-- '+name+'_END -->'
    left, rest = text.split(start, 1)
    _, right = rest.split(end, 1)
    return left+start+'\n'+content+'\n'+end+right


if __name__ == '__main__':
    environment_path = RESULTS/'environment.json'
    environment = read(environment_path)
    environment['packages'] = {name: version(name) for name in ('matplotlib', 'numpy', 'pytest', 'pillow')}
    environment['optimized_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                       for p in (ROOT/'rnai_designer').glob('*.py')}
    environment['benchmark_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                       for p in (ROOT/'benchmarks').glob('*.py')}
    environment_path.write_text(json.dumps(environment, indent=2), encoding='utf-8')
    comparisons = read(RESULTS/'comparison.json')
    assert len(comparisons) == 4
    lines = ['| nt | Ventanas | Total antes → después (s) | Factor | Cálculo antes → después (s) | Exportación antes → después (s) | Pico RSS antes → después (MiB) | CSV antes → después (MB) |',
             '|---:|---:|---:|---:|---:|---:|---:|---:|']
    for row in comparisons:
        a, b = row['before'], row['after']
        lines.append(f"| {row['nt']} | {a['windows']:,} | {a['total_s']:.2f} → {b['total_s']:.2f} | {row['speedup']:.2f}× | "
                     f"{a['compute_s']:.2f} → {b['compute_s']:.2f} | {a['export_s']:.2f} → {b['export_s']:.2f} | "
                     f"{a['peak_rss_bytes']/2**20:.1f} → {b['peak_rss_bytes']/2**20:.1f} | {a['csv_bytes']/1e6:.2f} → {b['csv_bytes']/1e6:.2f} |")
    largest = comparisons[-1]
    lines.append(f"\nEn 5000 nt, el pico RSS disminuyó {largest['peak_rss_reduction_percent']:.2f} %. "
                 'Exportación incluye CSV, JSON de selección y figuras. MiB = 2²⁰ bytes; MB = 10⁶ bytes.')
    observation = RESULTS/'reference_5000_memory_observation.json'
    if observation.is_file():
        memory = read(observation)
        lines.append(f"\nDurante la exportación de referencia de 5000 nt se observó además un máximo acumulado de "
                     f"{memory['PeakPagedMemorySize64']/1e9:.2f} GB de memoria comprometida, registrado en "
                     '`reference_5000_memory_observation.json`. Es un límite inferior del pico final, distinto del RSS. '
                     f"La ejecución optimizada final registró {largest['after']['peak_commit_bytes']/1e9:.2f} GB de compromiso máximo.")
    doc = ROOT/'docs'/'PERFORMANCE_OPTIMIZATION.md'
    text = replace_section(doc.read_text(encoding='utf-8'), 'BENCHMARK_TABLE', '\n'.join(lines))
    test_log = RESULTS/'tests.txt'
    raw_tests = test_log.read_bytes()
    tests = raw_tests.decode('utf-16' if raw_tests.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8-sig')
    test_log.write_text(tests, encoding='utf-8')
    summary = next(line for line in reversed(tests.splitlines()) if 'passed' in line)
    validation = [summary.strip('= '), '',
                  'La omisión corresponde exclusivamente al headless Conda legado. La prueba de Bowtie real usa los ejecutables empaquetados y un índice sintético temporal.', '']
    for mode in ('before', 'after'):
        trace = read(RESULTS/mode/'500'/'measurement_traced.json')
        validation.append(f"- `tracemalloc`, {mode}, 500 nt: pico Python {trace['tracemalloc_peak_bytes']/2**20:.2f} MiB; "
                          f"tiempo instrumentado {trace['total_s']:.2f} s (excluido de la tabla principal).")
    probe = read(RESULTS/'parallel_probe.json')
    validation.extend(['', f"Sonda Windows, dos referencias sintéticas de 500 nt: secuencial {probe['sequential_s']:.2f} s; "
                       f"spawn, cálculo y transferencia {probe['spawn_and_transfer_s']:.2f} s; factor {probe['speedup']:.2f}×. "
                       'Se verificaron todas las filas y el ganador. No se habilitó paralelismo del producto: esta sonda no valida el rendimiento ni el aislamiento de Conda/RNAplfold simultáneos.', '',
                       '| nt | Bowtie tres especies, primera búsqueda (s) | Repetición anterior (s) | Reutilización verificada (s) | RNAplfold directo (s) |',
                       '|---:|---:|---:|---:|---:|'])
    for row in read(RESULTS/'native_tools.json'):
        validation.append(f"| {row['nt']} | {row['rnai_designer_three_species_first_s']:.3f} | "
                          f"{row['reference_same_gene_second_reference_s']:.3f} | {row['rnai_designer_same_gene_second_reference_s']:.3f} | {row['rnaplfold_direct_s']:.3f} |")
    validation.append('\nLa primera búsqueda incluye fingerprints, versión, ejecución y parseo. El JSON conserva también el tiempo aislado de `execute` y el número de invocaciones: tres al buscar y cero al reutilizar. La caché puede añadir costo a la primera búsqueda; su beneficio está en las repeticiones. Las herramientas nativas no se aceleraron ni se modificaron: se evitó repetir búsquedas completas. Se verificó igualdad de los criterios científicos en ambos modos. Los índices son pequeños y sintéticos; estos tiempos no estiman búsquedas en las bases biológicas instaladas.')
    text = replace_section(text, 'VALIDATION', '\n'.join(validation))
    doc.write_text(text, encoding='utf-8')
    print('Updated docs/PERFORMANCE_OPTIMIZATION.md from four measured cases.')
