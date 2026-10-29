"""Configuración portable e inventario de índices; nunca escribe en las bases."""
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import os
import shutil
import subprocess

from .data_models import read_json

INDEX_PARTS = (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")
DEFAULT_DATABASES = {"Dalbulus maidis": "Dalbulus_maidis", "Homo sapiens": "Homo_sapiens",
                     "Zea mays": "zea_mays", "Apis mellifera": "Apis_mellifera"}


def configured_path(value: str, base: Path) -> Path:
    expanded = os.path.expandvars(value)
    path = Path(expanded).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


@dataclass
class Settings:
    database_dir: Path
    bowtie: Path | None
    bowtie_inspect: Path | None
    rnaplfold: Path | None
    databases: dict[str, dict]

    @classmethod
    def load(cls, config: Path | None = None, database_dir: Path | None = None):
        data = read_json(config) if config else {}
        if not isinstance(data, dict) or data.get("schema_version", "1.0") != "1.0":
            raise ValueError("Configuración: schema_version debe ser 1.0")
        base = config.resolve().parent if config else Path.cwd()
        local = os.environ.get("LOCALAPPDATA")
        home = os.environ.get("SIFI_HOME") or (str(Path(local) / "siFi2015") if local else None)
        db = database_dir or os.environ.get("SIFI_DB_DIR") or data.get("database_dir")
        if db is None:
            if home is None:
                raise ValueError("Definí --db-dir, SIFI_DB_DIR o LOCALAPPDATA")
            db = str(Path(home) / "Databases")
        db_path = configured_path(str(db), base)
        home_path = configured_path(home, base) if home else db_path.parent

        def tool(key, folder, executable):
            explicit = os.environ.get("SIFI_" + key.upper()) or data.get(key)
            if explicit:
                return configured_path(explicit, base)
            found = home_path / folder / (executable + (".exe" if os.name == "nt" else ""))
            return found if found.is_file() else (Path(shutil.which(executable)) if shutil.which(executable) else None)

        databases = {}
        configured = data.get("databases", {})
        if not isinstance(configured, dict):
            raise ValueError("databases debe ser objeto por organismo")
        for organism, default in DEFAULT_DATABASES.items():
            entry = configured.get(organism, {"prefix": default})
            if not isinstance(entry, dict) or not isinstance(entry.get("prefix"), str):
                raise ValueError(f"Base de {organism}: falta prefix")
            prefix = entry["prefix"]
            if not prefix or prefix in (".", "..") or "/" in prefix or "\\" in prefix:
                raise ValueError("prefix debe ser un nombre de índice dentro de database_dir")
            databases[organism] = dict(entry)
        return cls(db_path, tool("bowtie", "Bowtie", "bowtie"),
                   tool("bowtie_inspect", "Bowtie", "bowtie-inspect"),
                   tool("rnaplfold", "RNAplfold", "RNAplfold"), databases)


def index_files(directory: Path, prefix: str, extension: str = ".ebwt") -> list[Path]:
    return [directory / (prefix + part + extension) for part in INDEX_PARTS]


def inventory(settings: Settings, inspect: bool = False, timeout: float = 30) -> dict:
    directory = settings.database_dir
    result = {"schema_version": "1.0", "directory": str(directory), "exists": directory.is_dir(),
              "files": [], "indexes": [], "warnings": []}
    if not directory.is_dir():
        result["warnings"].append("Directorio de bases ausente; no se descargó ni creó")
        return result
    groups = set()
    # Los paquetes locales pueden guardar un segundo juego en Off_Targets/.
    # Inventariar un nivel adicional, sin mezclarlo con los prefijos activos.
    paths = sorted(p for p in directory.iterdir() if p.is_file())
    paths.extend(sorted(p for folder in directory.iterdir() if folder.is_dir()
                        for p in folder.iterdir() if p.is_file()))
    for path in paths:
        if not path.is_file():
            continue
        stat = path.stat()
        result["files"].append({"name": path.relative_to(directory).as_posix(), "size": stat.st_size,
                                "modified_at": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()})
        for extension in (".ebwt", ".ebwtl"):
            for part in sorted(INDEX_PARTS, key=len, reverse=True):
                suffix = part + extension
                if path.name.endswith(suffix):
                    groups.add((path.parent, path.name[:-len(suffix)], extension))
                    break
    for index_directory, prefix, extension in sorted(groups, key=lambda item: (str(item[0]), item[1], item[2])):
        parts = index_files(index_directory, prefix, extension)
        missing = [p.name for p in parts if not p.is_file() or p.stat().st_size == 0]
        entry = {"prefix": prefix, "format": "Bowtie 1" + (" large" if extension == ".ebwtl" else ""),
                 "directory": index_directory.relative_to(directory).as_posix(),
                 "extension": extension, "complete": not missing, "missing": missing,
                 "organism_labels": [name for name, db in settings.databases.items()
                                     if index_directory == directory and db["prefix"] == prefix],
                 "identity_note": "Asociación por configuración/nombre; no certifica taxonomía ni versión del ensamblado"}
        if inspect and not missing and settings.bowtie_inspect:
            try:
                process = subprocess.run([str(settings.bowtie_inspect), "-n", str(directory / prefix)],
                                         capture_output=True, text=True, errors="replace", timeout=timeout,
                                         creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
                names = process.stdout.splitlines()
                entry["inspection"] = {"returncode": process.returncode,
                                       "reference_count": len(names) if process.returncode == 0 else None,
                                       "sample_ids": names[:5], "stderr": process.stderr[:2000]}
            except (OSError, subprocess.TimeoutExpired) as error:
                entry["inspection"] = {"error": str(error)}
        result["indexes"].append(entry)
    return result


def doctor(settings: Settings) -> dict:
    tools = {}
    for name in ("bowtie", "bowtie_inspect", "rnaplfold"):
        path = getattr(settings, name)
        entry = {"path": str(path) if path else None, "available": False}
        if path:
            try:
                command = [str(path), "--help" if name == "bowtie_inspect" else "--version"]
                process = subprocess.run(command, capture_output=True, text=True, errors="replace", timeout=10,
                                         creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
                entry.update(available=process.returncode == 0, returncode=process.returncode,
                             version_output=(process.stdout + process.stderr)[:1500])
            except (OSError, subprocess.TimeoutExpired) as error:
                entry["error"] = str(error)
        tools[name] = entry
    from .legacy_runtime import detect_sifi_runtime
    sifi_root = Path(__file__).resolve().parents[2] / "siFi21--1.2.3"
    runtime = detect_sifi_runtime(sifi_root)
    return {"tools": tools, "databases": inventory(settings),
            "designer": {"automatic": bool(runtime.get("available")),
                         "headless_adapter_available": bool(runtime.get("available")),
                         "runtime": runtime,
                         "reason": ("Adaptador headless comprobado con los métodos originales de si-Fi. "
                                    "La selección de main target requiere --main-target o interacción en consola."
                                    if runtime.get("available") else runtime.get("reason"))}}
