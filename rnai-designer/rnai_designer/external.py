"""Procesos externos sin shell; conserva salida, errores y terminación."""
from datetime import datetime, timezone
from pathlib import Path
import subprocess


def execute(command: list[str], directory: Path, timeout: float) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    result = {"command": command, "started_at": datetime.now(timezone.utc).isoformat(),
              "status": "error", "returncode": None}
    with (directory / "stdout.log").open("wb") as stdout, (directory / "stderr.log").open("wb") as stderr:
        try:
            process = subprocess.run(command, cwd=directory, stdout=stdout, stderr=stderr, timeout=timeout,
                                     creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            result.update(returncode=process.returncode,
                          status="completed" if process.returncode == 0 else "error")
        except subprocess.TimeoutExpired:
            result.update(status="timeout", error=f"Tiempo máximo excedido: {timeout} s")
        except OSError as error:
            result["error"] = str(error)
    result["finished_at"] = datetime.now(timezone.utc).isoformat()
    return result
