# Instalación, herramientas y reproducibilidad

[Índice](README.md) · [Batch](batch-y-reportes.md) · [Créditos](creditos-y-licencias.md)

## Versiones verificadas

Se distingue requisito declarado, metadatos de esta máquina y evidencia histórica. No se ejecutaron herramientas científicas para obtener esta tabla.

| Componente | Versión/requisito | Evidencia |
|---|---|---|
| rnai-designer | **0.3.0** | [`pyproject.toml`](../rnai-designer/pyproject.toml) y [`__init__.py`](../rnai-designer/rnai_designer/__init__.py) |
| Python de rnai-designer | Requiere **>=3.10**; inspeccionado **3.14.10** | Declaración y lectura de `sys.version`; sin matriz completa de compatibilidad |
| Matplotlib Python 3 | Rango **>=3.8,<4**; instalado **3.11.2** | Dependencia obligatoria y metadatos locales |
| NumPy / Pillow Python 3 | **2.5.4 / 13.0.0** instalados | Metadatos; no son dependencias directas declaradas del paquete |
| pytest | Extra **>=8,<10**; instalado **9.2.0** | Opcional para pruebas |
| si-Fi original | **v1.2.3-0008** | [README original](../siFi21--1.2.3/README.md) |
| Python legado | **2.7.13**, **64 bits** en el manifiesto observado | `conda-meta` y `run.json` histórico de BicC_Dm4 |
| Biopython / NumPy legados | **1.76 / 1.16.5** | Metadatos del entorno `sifi_32bit` |
| Matplotlib / SciPy / Seaborn legados | **2.2.5 / 1.2.3 / 0.9.1** | Metadatos de distribución, sin importar módulos |
| PyQt / Qt / SIP legados | **4.11.4 / 4.8.7 / 4.18** | Archivos `conda-meta` |
| Bowtie | **0.11.3** registrado | Manifiesto histórico `offtarget_runs/human/run.json` de BicC_Dm4; familia clásica, no Bowtie 2 |
| RNAplfold / ViennaRNA | **Versión exacta pendiente** | Ejecutable presente; metadatos Windows de versión vacíos. No se ejecutó `--version` |
| Conda | **Versión exacta pendiente** | No se invocó Conda durante la revisión |
| Bases biológicas | **Versiones/ensamblados pendientes** | Configuración de ejemplo con `source`/`version` nulos; manifiesto Bowtie consultado con `unknown` |

`sifi_32bit` es un nombre y no certifica arquitectura de 32 bits. Los metadatos instalados describen esta máquina; no prueban que todas las importaciones o herramientas funcionen correctamente. El README original identifica la versión de si-Fi; el nombre de la carpeta no basta.

Fuentes de dependencias Python 3: [requisitos básicos](../rnai-designer/requirements.txt), [desarrollo](../rnai-designer/requirements-dev.txt), [versiones fijadas](../rnai-designer/requirements-tested.txt). No se encontró un lockfile completo y verificado del entorno legado.

## Python 3 en Windows/PowerShell

Desde la raíz Git:

```powershell
cd rnai-designer
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
```

No hace falta activar el entorno ni modificar la política de PowerShell. Opcionalmente, para desarrollo:

```powershell
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
```

Para reproducir las dependencias fijadas de la entrega, utilizando un Python compatible:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-tested.txt
```

Los rangos no garantizan ruedas binarias para todas las versiones de Python. No se verificó una matriz completa de plataformas.

## Entorno legado y recursos externos

El motor requiere Conda y Python 2.7 con PyQt4, Biopython que conserve `Bio.Alphabet`, NumPy, Matplotlib, SciPy y Seaborn compatibles. Sus versiones modernas no son sustitutos directos. `siFi21--1.2.3/setup.py` es una receta histórica de empaquetado con py2exe, no un instalador reproducible del entorno headless.

**Pendiente:** publicar una exportación verificable del entorno legado, con canales, builds y procedencia de binarios. No se proporciona un supuesto comando `conda create`: falta una receta completa validada.

Se requieren Bowtie clásico, `bowtie-inspect`, RNAplfold y bases para *D. maidis*, *H. sapiens*, *Z. mays* y *A. mellifera*. `bowtie-build` es necesario si se construyen índices o se ejecutan pruebas de integración que lo utilicen. Ni estas herramientas ni las bases se descargan automáticamente al instalar el paquete Python.

## Configuración y rutas

Valores predeterminados reales del programa en Windows:

```text
%LOCALAPPDATA%\siFi2015\Bowtie\bowtie.exe
%LOCALAPPDATA%\siFi2015\Bowtie\bowtie-inspect.exe
%LOCALAPPDATA%\siFi2015\RNAplfold\RNAplfold.exe
%LOCALAPPDATA%\siFi2015\Databases\
```

El código original se reconoce por `sifi_pipeline.py`, `general_helpers.py` y `show_plot.py`. Conda se busca en PATH y ubicaciones conocidas; se prioriza `sifi_32bit` entre los entornos cuyo nombre contiene `sifi`.

| Variable | Función |
|---|---|
| `SIFI_ROOT` | Directorio de fuentes originales |
| `SIFI_CONDA_ENV` | Nombre o prefijo del entorno Python 2.7 |
| `SIFI_HOME` | Raíz de recursos |
| `SIFI_DB_DIR` | Directorio de índices |
| `SIFI_BOWTIE`, `SIFI_BOWTIE_INSPECT`, `SIFI_RNAPLFOLD` | Ejecutables, no carpetas |
| `CONDA_EXE` | Localización explícita de Conda |

Desde `rnai-designer/`, para la disposición de este repositorio:

```powershell
$env:SIFI_ROOT = (Resolve-Path ..\siFi21--1.2.3).Path
$env:SIFI_CONDA_ENV = "sifi_32bit"
```

El entorno del segundo comando debe existir; en otra máquina puede tener otro nombre. No ejecutar instalaciones ni modificar el entorno de un análisis activo.

La [configuración de ejemplo](../rnai-designer/examples/config.windows.json) utiliza los prefijos `Dalbulus_maidis`, `Homo_sapiens`, `zea_mays` y `Apis_mellifera`. Las rutas del JSON son relativas a ese JSON; las rutas CLI son relativas al directorio de trabajo. Para los ajustes correspondientes, prevalecen las opciones CLI sobre variables y estas sobre JSON.

Se esperan índices activos directamente en `database_dir`; el inventario de una subcarpeta no la activa. Se reconocen seis partes de índices `.ebwt`/`.ebwtl`; revise el diagnóstico y la compatibilidad efectiva del ejecutable para cada formato.

## Doctor y diagnóstico

Ejecutar antes de iniciar o después de terminar el batch activo:

```powershell
.\.venv\Scripts\python.exe -m rnai_designer doctor
.\.venv\Scripts\python.exe -m rnai_designer doctor --config examples/config.windows.json
.\.venv\Scripts\python.exe -m rnai_designer databases --config examples/config.windows.json
```

`doctor` ejecuta comprobaciones externas: consulta herramientas, Conda e importaciones legadas, además de inventariar índices. No analiza un gen, pero no es una mera lectura de metadatos. Revise `available`, `reason`, errores y versiones en el JSON; el código de salida por sí solo no certifica disponibilidad total.

Para guardar un diagnóstico posterior:

```powershell
.\.venv\Scripts\python.exe -m rnai_designer doctor --out diagnostico.json
```

Elija un archivo nuevo si necesita conservar diagnósticos anteriores: `--out` puede escribir sobre el archivo indicado. Ninguno de estos comandos se ejecutó durante esta edición documental.

## Registro mínimo de reproducibilidad

Conservar commit y cambios locales, intérpretes y dependencias; FASTA y hash de entrada; procedencia de genes; parámetros; herramientas y builds; origen, versión y hashes de índices; JSON si-Fi, máscara, manifiestos, logs e intentos off-target.

Distinguir fixtures ficticios de `examples/`, pruebas sintéticas y análisis biológicos reales. Los reportes no reemplazan un registro completo del entorno. Una versión desconocida no debe inferirse del nombre, la fecha o el identificador GITV.
