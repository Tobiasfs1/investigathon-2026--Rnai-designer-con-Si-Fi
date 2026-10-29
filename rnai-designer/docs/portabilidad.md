# Instalación portable en Windows

`rnai-designer` usa Python 3.10 o posterior. El motor original si-Fi sigue
ejecutándose en Python 2.7 dentro de Conda. No se cambian sus algoritmos ni bases.

## Preparar otra computadora

1. Copiar `rnai-designer` y el código original si-Fi como carpetas hermanas.
   El nombre de la carpeta original puede cambiar: se reconoce por
   `sifi_pipeline.py`, `general_helpers.py` y `show_plot.py`.
2. Instalar Python 3 y crear un entorno nuevo; no copiar `.venv` de otra PC:

   ```powershell
   python -m venv .venv
   .\.venv\Scripts\python.exe -m pip install -e ".[test]"
   ```

3. Instalar Miniconda, Miniforge o Anaconda y restaurar el entorno antiguo
   `sifi_32bit`, o configurar otro entorno equivalente. Debe contener Python 2.7,
   PyQt4, Biopython con `Bio.Alphabet`, NumPy, Matplotlib, Seaborn y SciPy
   compatibles con Python 2. Instalar las versiones actuales de esos paquetes
   no reemplaza al entorno antiguo. `doctor` muestra las versiones detectadas.
   El nombre `sifi_32bit` no impone arquitectura: se informa la arquitectura real.
4. Instalar Bowtie **1**, incluyendo `bowtie-inspect` y, para las pruebas de
   integración, `bowtie-build`; instalar RNAplfold de ViennaRNA. Por defecto se
   buscan `Bowtie/bowtie.exe`, `Bowtie/bowtie-inspect.exe` y
   `RNAplfold/RNAplfold.exe` dentro de `%LOCALAPPDATA%/siFi2015`.
5. Copiar/instalar los índices Bowtie 1 en `%LOCALAPPDATA%/siFi2015/Databases`.
   Cada índice debe tener las seis partes no vacías. Los prefijos por organismo
   se pueden ajustar en `databases` del JSON de configuración existente.

## Detección y configuración

Conda se busca primero en PATH (`conda.exe`, luego `conda`), después en
`CONDA_EXE`, en las raíces deducibles de `CONDA_PREFIX` y en instalaciones
habituales bajo el perfil del usuario, LOCALAPPDATA, PROGRAMDATA y ProgramFiles.
Incluye Miniconda, Miniforge, Anaconda y Mambaforge; no depende de un usuario ni
de una unidad concretos. Para una instalación arbitraria fuera de esas raíces,
poner Conda en PATH o definir `CONDA_EXE`.

Se consulta `conda env list --json` y se prioriza `sifi_32bit` entre los entornos
con `sifi` en el nombre. Una selección explícita admite nombre o prefijo. Cada
candidato se ejecuta mediante `conda run -p`: se comprueba Python 2.7 y cada
dependencia e importación del puente. No se importa `main`, que corresponde a
la GUI y puede exigir componentes ajenos a la ejecución headless.

Las opciones de línea de comandos prevalecen sobre las variables de entorno;
estas prevalecen sobre los campos JSON. Las rutas del JSON son relativas al
archivo de configuración. Las variables que pasan por `Settings` también usan
esa base cuando se proporciona `--config`; sin configuración usan el directorio
de trabajo. Las rutas de CLI son relativas al directorio de trabajo.
Se expanden variables de entorno y `~`.

| Componente | Variable | Campo JSON | Opción CLI |
|---|---|---|---|
| Código original | `SIFI_ROOT` | `sifi_root` | `--sifi-root` |
| Entorno Python 2 | `SIFI_CONDA_ENV` | `legacy_conda_env` | `--legacy-conda-env` |
| Directorio común de recursos | `SIFI_HOME` | `sifi_home` | — |
| Índices | `SIFI_DB_DIR` | `database_dir` | `--db-dir` |
| Bowtie | `SIFI_BOWTIE` | `bowtie` | — |
| bowtie-inspect | `SIFI_BOWTIE_INSPECT` | `bowtie_inspect` | — |
| RNAplfold | `SIFI_RNAPLFOLD` | `rnaplfold` | — |

Las tres variables de herramientas indican **ejecutables**, no carpetas.
Sin configuración explícita, las herramientas también se buscan en PATH.
Si falta LOCALAPPDATA se usa `~/AppData/Local/siFi2015`.
El código original se busca desde la ubicación del paquete, independientemente
del directorio de trabajo, en la carpeta del proyecto, su padre y sus hermanas.
Si está fuera de esa estructura, configurar `SIFI_ROOT` o `sifi_root`.

Ejemplo de archivo portable ubicado en la raíz de `rnai-designer`:

```json
{
  "schema_version": "1.0",
  "sifi_root": "../original-sifi",
  "legacy_conda_env": "sifi_32bit",
  "sifi_home": "../recursos-sifi"
}
```

## Diagnóstico

Desde `rnai-designer`, con Python 3:

```powershell
.\.venv\Scripts\python.exe -m rnai_designer doctor
.\.venv\Scripts\python.exe -m rnai_designer doctor --config config.json
.\.venv\Scripts\python.exe -m rnai_designer doctor --legacy-conda-env otro_entorno --sifi-root ..\original-sifi
.\.venv\Scripts\python.exe -m rnai_designer doctor --out output\diagnostico.json
```

El JSON incluye Conda y su versión, herramientas y versiones, inventario de bases,
estado de cada índice configurado, prefijo e intérprete del entorno, versión y
arquitectura de Python, versiones de dependencias, rutas de módulos y errores
individuales. Comprobar `designer.runtime.available`, los `tools.*.available` y
`databases.configured_indexes[*].available`. Que el motor esté disponible no
significa que las herramientas externas estén instaladas.

Por compatibilidad, `doctor` devuelve 0 cuando pudo generar el diagnóstico,
incluso si informa componentes ausentes; errores de configuración devuelven 2.
El diagnóstico no instala ni modifica Conda, el código si-Fi o las bases.

El mensaje anterior mezclaba descubrimiento, dependencias e importaciones en un
solo error: encontrar un entorno e importar dos módulos manualmente no probaba
que la validación completa del programa fuera correcta. Ahora se informa la
causa y el candidato exactos; no se atribuye todo fallo a la ausencia de Conda.
En la computadora revisada se reprodujo el fallo original: `main` importa
`imageviewer`, que requiere `PIL`, ausente en ese entorno. Los módulos del
análisis sí importan correctamente. La GUI continúa requiriendo sus propias
dependencias; la disponibilidad headless se comprueba por separado.

## Pruebas

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

Si Windows deniega acceso al directorio temporal global de pytest, usar una
carpeta temporal nueva y exclusiva dentro del proyecto:

```powershell
.\.venv\Scripts\python.exe -m pytest -q --basetemp=.test-portability-local
```

Las pruebas de integración necesitan las herramientas externas instaladas;
se omiten cuando no están disponibles. Los tests de portabilidad simulan
usuarios, ubicaciones, configuración y fallos de dependencias sin modificar
instalaciones reales.
