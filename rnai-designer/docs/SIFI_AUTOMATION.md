# Automatización del flujo original de si-Fi

## Cómo se calculan y exponen las regiones verdes

La inspección corresponde a la copia local `siFi21--1.2.3` (README: v1.2.3-0008). Los límites verdes exactos no están serializados como `green_regions` en el JSON temporal, pero la máscara sí existe explícitamente durante la ejecución:

- `sifi_pipeline.py:SifiPipeline.create_sirnas` (líneas 286–304) crea cada ventana y la nombra `sirna<start+1>`. `data_to_json` (líneas 184–263) extrae el sufijo como `sirna_position`; es el inicio 1-based de la query.
- `general_helpers.py:get_target_data` (líneas 125–157) crea `plot_range = range(sirna_position, sirna_position + sirna_size)`. Devuelve `main_target_dict`, un diccionario `hit_name -> set(posiciones)`, además de las estructuras de off-target y eficiencia. La rama `if data['is_off_target']` usa truthiness de Python: `True` suma posiciones rojas; `False` y `null` entran en la rama main-target. `null` debe advertirse aunque se reproduzca ese comportamiento.
- `show_plot.py:DrawPlot.plot_design` (líneas 108–145) obtiene el diccionario del helper, une sus valores en `main_targets_plot` y dibuja un rectángulo verde por posición. El diccionario, por tanto, permite recuperar la máscara exacta antes de interpretar píxeles.
- `SifiPipeline.run_pipeline` (líneas 143–168) genera primero el JSON temporal, invoca el gráfico y devuelve `main_target_dict` como cuarto elemento de su tupla. La estructura está disponible en memoria al terminar la llamada. `main.py:start_sifi_pipeline` (líneas 144–153) recibe esa tupla, pero la GUI no exporta explícitamente el diccionario.
- `create_plots.py:create_design_plot` (líneas 44–60) repite el mismo acceso al helper en la ruta alternativa de gráficos. También es posible extraer la máscara llamando directamente a `get_target_data` una vez generado el JSON, antes de dibujar.

La unión puede tener varios runs separados: la GUI dibuja cada posición individual y los gaps quedan vacíos. Las líneas de texto `REGION VERDE` en `show_plot.py` (177–195) y `create_plots.py` (113–128), y el script `region_verde.py` (26–48), resumen solo `min(mask)..max(mask)`. Ese intervalo envolvente puede abarcar gaps y no es una enumeración exacta de regiones. Esas rutinas también tratan como 0-based valores derivados de `sirna_position`, aunque el pipeline los genera como 1-based, y desplazan el corte de secuencia. No usar esos textos/FASTA como fuente canónica de coordenadas.

El JSON temporal de `SifiPipeline.run_pipeline` es una lista de registros siRNA y se crea con `tempfile.mkstemp()` en el temporal predeterminado de Python; la ruta se devuelve en la tupla, pero no hay exportación persistente desde la GUI. El módulo `rnai-designer` acepta un sidecar lossless `green_mask.json` con `positions_1based`, exportado desde el diccionario original, y comprueba que coincida con la reconstrucción de filas. El puente `legacy_bridge/export_green_mask.py` llama al `general_helpers.get_target_data` original sin renderizar la imagen; debe ejecutarse con el Python 2.7 y Biopython antiguos que requiere si-Fi. La alternativa en Python 3 reproduce la misma semántica desde la lista JSON. El JSON temporal no contiene por sí mismo las regiones ni `main_target_dict` como campos explícitos.

## Entorno legado y ejecucion headless comprobada

En esta computadora, `conda env list` encontro las instalaciones base `C:\ProgramData\miniforge3` y `C:\Users\lucas\miniconda3`, ademas de `C:\Users\lucas\.conda\envs\sifi_32bit`. El nombre de este ultimo no coincide con su arquitectura: la sonda encontro Python 2.7.13 de **64 bits**, PyQt4 4.11.4, Qt 4.8.7, Biopython 1.76, NumPy 1.16.5, Matplotlib 2.2.5, Seaborn 0.9.1 y SciPy 1.2.3. `main`, `sifi_pipeline`, `general_helpers`, `show_plot` y `free_energy` importan correctamente desde el codigo original.

La GUI original puede iniciarse con `conda run -n sifi_32bit python "..\siFi21--1.2.3\main.py"`. Bowtie 0.11.3 y RNAplfold 2.1.5 estan instalados junto a si-Fi en `%LOCALAPPDATA%\siFi2015`, fuera del entorno Conda. El adaptador les pasa sus rutas absolutas. No se instalaron paquetes ni se modifico si-Fi, sus indices o bases.

La GUI no es necesaria para los calculos: `SifiPipeline.create_sirnas`, `run_rnaplfold`, `bowtie_to_lst` y `data_to_json` se pueden llamar sin instanciar ventanas. El impedimento era el selector PyQt4 de `get_main_target()` y el `plot.show()` posterior. `legacy_bridge/run_sifi_headless.py` invoca metodos originales y `general_helpers.get_target_data`, omitiendo ambas acciones de interfaz y exportando la mascara exacta. Se ejecuto exitosamente bajo el entorno encontrado contra un indice Bowtie **sintetico temporal**: calculo 160 registros de una secuencia de 180 nt, con metricas originales y mascara 1-180. Esto demuestra ejecucion headless para esta copia/instalacion; no valida eficacia biologica.

Con `workflow --run-sifi-headless` sin `--main-target`, primero se enumeran los IDs alineados por Bowtie y se ejecuta una corrida independiente de si-Fi para cada ID. `--main-target ID` restringe a una referencia; se puede repetir. Los archivos JSON y la mascara verde de cada target quedan separados. No se combinan coordenadas entre referencias.

Ejecucion automatica de todos los main targets detectados desde una sola secuencia:

```powershell
python -m rnai_designer workflow --gene-id ATPase --sequence-file inputs/ATPase.fasta --out output/ATPase_multi --run-sifi-headless
```

Para listar los IDs disponibles antes del analisis, agregar `--list-targets-only`. Para analizar una referencia especifica, agregar `--main-target GITV01001458_1`; para seleccionar varias, repetir esa opcion. Sin ella se procesan todas. Cada resultado queda en `references/NN_ID/`; `reference_comparison.json` y `all_references_windows.csv` consolidan las referencias sin mezclar sus coordenadas. Los parametros originales se guardan por referencia en `designer_handoff/run.json`.

## Una sola corrida para las etapas disponibles

`workflow` acepta el gen y su secuencia y con `--run-sifi-headless` enumera todos los main targets, ejecuta el adaptador original una vez por referencia y conserva JSON y mascara verdes separados. Una referencia fallida queda registrada en `reference_error.json` y las demas continuan. El reporte de raiz compara los runs verdes y ventanas en `all_references_windows.csv`, con una columna `main_target` para conservar la referencia de cada coordenada.

Preparación seguida de reanudación manual:

```powershell
python -m rnai_designer workflow --gene-id DM_gen --sequence-file gen.fasta --reference-source "repositorio" --reference-version "accesion" --out output/DM_gen
python -m rnai_designer workflow --resume output/DM_gen --sifi auto --config examples/config.windows.json
```

Una corrida que queda abierta esperando mientras se usa la GUI:

```powershell
python -m rnai_designer workflow --gene-id DM_gen --sequence-file gen.fasta --out output/DM_gen --wait-for-sifi --config examples/config.windows.json
```

Para extraer la máscara con el helper original durante esa corrida, agregar `--legacy-python C:/ruta/python2.exe`. También puede exportarse por separado desde un JSON existente:

```powershell
python C:/ruta/python2.exe rnai_designer/legacy_bridge/export_green_mask.py --sifi-json salida_sifi.json --sirna-size 21 --sifi-root ../siFi21--1.2.3 --out green_mask.json
```

Después, pasar `--green-mask green_mask.json` a `workflow` para conservar ese objeto y su procedencia en el reporte.

Antes de la corrida real, completar `designer_handoff/run.json` con el query ID y los parámetros/base/versión efectivamente usados. Para que la recuperación automática sea inequívoca, usar el FASTA de handoff y mantener el query ID del manifiesto. Si hay varios JSON compatibles, indicar la ruta explícita con `--sifi`.

El resultado es `pipeline_result.json`, junto con `integrated/report.json`, `offtarget_runs/` y `window_selection/`. El filtro se reporta por especie como `criterion_met`, `exact_match_detected`, `indeterminate`, `error` o `not_evaluated`. `criterion_met` significa que los datos importados acreditan la ausencia de coincidencias exactas de 21 nt bajo el alcance y parámetros documentados; no significa seguridad. La ventana ganadora conserva el puntaje de **calidad de evidencia**, no una puntuación de eficacia. GC, secuencias y etiquetas de si-Fi se presentan en campos independientes.

En `show_plot.py:141-145`, cada valor de la máscara se usa directamente como x y se dibuja el rectángulo `[p,p+1)` en un eje de `0` a longitud de query; no se transforma el inicio de JSON a `p-1`. Por eso las coordenadas de máscara/JSON y los cortes 0-based de la secuencia deben registrarse separadamente. `region_verde.txt` añade otra conversión incorrecta al sumar 1 a `min` y `max`.


### Multiples main targets y fidelidad de la mascara

En modo headless, Bowtie primero enumera todos los `hit_name`; si no se pasa `--main-target`, el workflow corre si-Fi una vez por target y conserva salidas en `references/01_ID/`, `references/02_ID/`, etc. `reference_comparison.json` incluye estado/error y resumen por referencia; `all_references_windows.csv` agrega `main_target` en cada fila. Un error de una referencia se escribe en su carpeta como `reference_error.json` y no detiene las demas.

Cada adaptador consulta `general_helpers.get_target_data` en su JSON original, por separado. Si una posicion siRNA aparece tanto en el target seleccionado como en alineamientos off-target, la implementacion original resta posiciones off-target al acumulado verde; en secuencias muy parecidas eso puede acortar o suprimir runs. El sidecar registra `main_offtarget_overlap_start_count` y una advertencia para que el resultado pueda auditarse. El workflow no cambia esa regla ni combina coordenadas entre IDs.
