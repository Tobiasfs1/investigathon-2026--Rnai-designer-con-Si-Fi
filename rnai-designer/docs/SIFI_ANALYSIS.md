# Análisis de la copia local de si-Fi

Inspección estática del 27/10/2026 de la carpeta hermana `siFi21--1.2.3/`, cuyo nombre local difiere de `si-Fi/`. No se modificó ni ejecutó RNAi Designer original. No se encontraron archivos de resultados reales dentro del código fuente. Las bases instaladas por separado se inspeccionaron después y están descritas en [SIFI_INTEGRATION.md](SIFI_INTEGRATION.md). El README declara **v1.2.3-0008** y licencia CC BY-NC-SA 2.0. El módulo nuevo no copia los algoritmos de predicción.

## Entorno y ejecución

El código requiere **Python 2**: `sifi_pipeline.py:SifiPipeline.run_pipeline` contiene `print x[0]` y `general_helpers.py:create_gbk` usa `iteritems`. `region_verde.py` también se declara Python 2. Los nueve `.pyc` del nivel principal tienen cabecera `03 f3 0d 0a`, correspondiente al magic 62211 de **CPython 2.7**, contrastado con [Python/import.c de CPython 2.7](https://github.com/python/cpython/blob/2.7/Python/import.c). Esto confirma que los bytecodes fueron compilados para esa familia; no revela la revisión 2.7.x ni garantiza que reflejen los fuentes actuales. No hay versiones fijadas de dependencias.

`main.py:main` crea QApplication y abre una ventana. Ejecución prevista desde la carpeta original: `python main.py`, usando un intérprete Python 2 compatible y dependencias legadas. El código importa PyQt4 (QtGui, QtCore, QtNetwork), FileDialog, Biopython, NumPy, matplotlib, seaborn y Pillow. `Bio.Alphabet` y `np.float` obligan a revisar versiones antiguas. `setup.py` usa py2exe y pytz e incluye sip, componentes SciPy, zmq y backend Tk para empaquetado: esto no significa que cada componente participe en la predicción.

`ToCopy/` contiene ejecutables Windows de Bowtie y RNAplfold con DLL. `main.py:MyMainWindow.__init__` prepara carpetas bajo datos locales `siFi2015/`; `general_helpers.copying_files` copia recursos desde `ToCopy/` relativo al directorio de ejecución. `database_helpers.create_bowtie_database` llama `bowtie-build` para indexar FASTA. No se reconstruyó el entorno legado ni se certificaron versiones de esos ejecutables.

## Flujo y eficiencia predicha

En `sifi_pipeline.py`, `SifiPipeline.create_sirnas` enumera ventanas solapadas del tamaño elegido, con nombres `sirna1...` e inicio 1-based. `run_bowtie` llama `bowtie -a -v <mismatches> -y <db> -f <fasta> <salida>`. El código no comprueba correctamente el retorno del proceso.

`run_rnaplfold` llama `RNAplfold -W 80 -L 40 -u <sirna_size> -T 22`, crea un directorio temporal, carga `<query_name>_lunp` y elimina las primeras `sirna_size-1` filas. `data_to_json` toma la fila `sirna_position-1` y columna `accessibility_window`. Se conserva ese valor como accesibilidad importada; no se interpreta como probabilidad de eficacia biológica.

`free_energy3` y `free_energy_dangling_ends` llaman `free_energy.calculate_free_energy`, que contiene tablas termodinámicas. `strand_selection` compara energías de extremos; `end_stability` compara su diferencia contra un umbral; `pair_probability` compara accesibilidad contra otro umbral. `check_efficient` combina los filtros habilitados. `calculate_efficiency` agrega reglas terminales que pueden alterar esa combinación. Por eso `is_efficient` no siempre equivale al AND de los tres booleanos. Los índices terminales utilizados son `[1]` y `[sirna_size-3]`, aunque algunos comentarios sobre posiciones resultan confusos. No se reimplementan esas reglas.

`main.py:test` deshabilita controles si se permiten mismatches. `start_sifi_pipeline` establece entonces `no_efficience`. En esa rama, `calculate_efficiency` devuelve `False` y los tres filtros `null`: ese False puede significar cálculo desactivado.

## JSON temporal y coordenadas verificadas

`run_pipeline` usa `tempfile.mkstemp()` sin `dir` y escribe `json.dump(json_lst, indent=4)` en `<ruta temporal>.json`. Es **una lista JSON**, no JSONL ni un objeto con metadatos. Está en el temporal predeterminado de Python, habitualmente TEMP en Windows; no necesariamente en `self.temp_location`. No se encontró borrado explícito del JSON en este flujo, pero el sistema puede limpiar temporales. Su nombre aleatorio no identifica gen ni base. La ruta se devuelve como segundo elemento de la tupla de `run_pipeline`. No hay exportación explícita del JSON en la interfaz inspeccionada.

Campos construidos en `SifiPipeline.data_to_json`:

| Campo | Contenido |
|---|---|
| `query_name`, `sirna_name` | ID FASTA y nombre de ventana |
| `sirna_position` | Inicio en la query, 1-based |
| `sirna_sequence` | Secuencia exportada por Bowtie |
| `is_efficient` | Clasificación del software, no validación experimental |
| `strand_selection`, `end_stability`, `target_site_accessibility` | Booleanos de filtros o null |
| `accessibility_value` | Valor numérico extraído de RNAplfold |
| `is_off_target` | Diseño: según selección manual; modo off-target: null |
| `hit_name` | ID de referencia alineada |
| `reference_strand_pos` | Offset Bowtie, 0-based |
| `strand` | Hebra del hit |
| `mismatches` | Descriptor textual Bowtie o entero 0, no necesariamente conteo |

No se exportan energías numéricas, `thermo_effcicient`, parámetros, versión, base, organismo ni evidencia experimental. `bowtie_to_lst` agrega entero 0 cuando hay siete columnas. Si hay descriptor, se conserva como texto.

Varias filas pueden corresponder al mismo siRNA y distintos hits. Para gráficos se dispone de inicio, longitud de secuencia, accesibilidad y etiquetas. El intervalo correcto sobre query es 1-based inclusivo `[inicio, inicio+longitud-1]`. El prototipo grafica valores por inicio de ventana, sin interpolar observaciones faltantes ni optimizar regiones.

## Modos, resultados y limitaciones observadas

`mode=0` llama `get_main_target`, que abre un popup para seleccionar los objetivos principales. Las otras referencias se marcan off-target. `mode=1` no solicita esa selección y exporta `is_off_target=null`.

`main.py:set_design_settings` establece tamaño 21, mismatches 0, filtros activados, umbrales 1 y 0.1, ventana 8. `set_offtarget_settings` desactiva end/accessibility y mantiene strand/terminal. Son defaults de la interfaz, no parámetros demostrados para cualquier archivo. `run_pipeline` ejecuta RNAplfold en ambos modos.

**Limitación crítica:** en `data_to_json`, tanto la construcción como el agregado de registros están dentro de `if strand == '+'`. Se omiten hits de hebra `-`. En la rama sin objetivos (`strand=None`) tampoco se agregan registros, pese a los comentarios que prometen graficar eficiencia sin hits. Off-target sin hits no genera JSON útil. Un archivo ausente o vacío nunca demuestra ausencia de riesgo. Además, `run_pipeline` retorna dentro del bucle FASTA: no se debe asumir procesamiento completo de múltiples queries.

`show_plot.py:DrawPlot.plot_design/plot_offtarget` genera PNG temporales. `create_plots.py:CreatePlots` contiene otro trazado, pero su importación en el pipeline está comentada. `show_plot.DrawPlot.export_table` e `imageviewer.ImageViewer.export_table` exportan CSV separado por `;` con `Targets;Total siRNA hits;Efficient siRNA hits`. Ese resumen pierde posiciones y métricas individuales. `imageviewer.export_gbk` y `general_helpers.create_gbk` permiten GenBank, pero el bloque que abre ImageViewer en `main.py:start_sifi_pipeline` está comentado: no se garantiza acceso a esa opción en la UI activa.

`general_helpers.get_target_data` acumula conjuntos entre targets; sus agrupaciones pueden mezclar referencias y depender del orden. `region_verde.py` incluye una CLI local en comentarios para extraer una región de un JSON existente, pero interpreta posiciones como 0-based donde el productor usa 1-based y toma min/max aunque haya discontinuidades. El prototipo no utiliza ese script.

## Ejecución externa y documentación disponible

No se encontró CLI general documentada de predicción ni API pública estable. El README local describe el programa de escritorio y su licencia, sin guía de instalación detallada. `main.py:show_help` enlaza un PDF remoto, no incluido ni usado como evidencia en esta auditoría. `region_verde.py` procesa resultados existentes; no es una CLI general de si-Fi.

`SifiPipeline` importa popup y show_plot; `run_pipeline` es una propiedad que abre ventanas y requiere interacción. Un adaptador Python 2 headless de RNAi Designer no se incluyó porque faltan una interfaz estable y selección automática validada de main targets. Python 3 ya prepara el FASTA y recupera/importa JSON de esa corrida manual. Para búsquedas exactas de off-target se añadió una integración separada con Bowtie local: [SIFI_INTEGRATION.md](SIFI_INTEGRATION.md). No se importa ningún módulo Python 2.
