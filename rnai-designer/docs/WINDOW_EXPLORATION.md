# Exploración de ventanas en regiones verdes de RNAi Designer

## Qué representa el verde en si-Fi

El gráfico activo se implementa en `siFi21--1.2.3/show_plot.py`, método `DrawPlot.plot_design`. El verde no indica eficiencia. El título de la leyenda `Main target` corresponde a posiciones cubiertas por hits del objetivo principal. Las listas de posiciones se construyen en `siFi21--1.2.3/general_helpers.py`, función `get_target_data`, a partir de `sirna_position`, `sirna_size`, `hit_name` e `is_off_target`.

Para cada fila, el helper crea `range(sirna_position, sirna_position + sirna_size)`: `sirna_position` se genera como inicio 1-based en `sifi_pipeline.py:create_sirnas`, así que el intervalo de bases es inclusivo `[inicio, inicio+size-1]`. Si `is_off_target` es truthy, agrega ese tramo a `off_target_positions`. Si es falsey, acumula el tramo en `main_target_positions` y asigna a `main_target_dict[hit_name]` una instantánea de `main_target_positions - off_target_positions`. **El `null` también es falsey en el código original y, por lo tanto, entra por la rama main-target**, aunque su significado sea ambiguo; se conserva una advertencia. El gráfico toma la unión de los valores de `main_target_dict` como `main_targets_plot` y pinta un rectángulo verde, ancho 1, por cada posición.

La máscara de posiciones sí existe explícitamente en memoria antes de renderizar: `general_helpers.get_target_data` devuelve `main_target_dict`, un diccionario `hit_name -> set(posiciones)`. `SifiPipeline.run_pipeline` pasa esa estructura a `DrawPlot`, la devuelve en su tupla y `show_plot.py:plot_design` la une para construir `main_targets_plot` antes de crear los parches. Por lo tanto se puede extraer directamente desde el pipeline o invocando el helper al leer el JSON, sin interpretar píxeles. El JSON temporal de si-Fi no serializa ese diccionario ni runs ya agrupados; el explorador acepta el sidecar lossless del helper original y, si no existe, reconstruye la misma unión desde las filas.

Tampoco es equivalente a `is_efficient`: si-Fi dibuja por separado una curva roja de siRNAs clasificados como eficientes y una curva azul de todos los hits del main target. El verde es la cobertura main-target. Una fila main-target puede tener `is_efficient == false`.

La implementación antigua guarda instantáneas conforme recorre las filas. `main_target_positions` y `off_target_positions` se acumulan; un off-target posterior no modifica instantáneas anteriores. La agrupación puede, por tanto, depender del orden del JSON y mezclar referencias. El explorador conserva el orden original y reproduce esa lógica, incluyendo el tratamiento falsey de `null`, y emite una advertencia. Para el resultado más fiel se necesita un JSON de **RNAi design**, no el modo Off-target Predictions.

`DrawPlot.plot_design` dibuja las marcas verdes primero y rectángulos rojos de off-target después. En el gráfico original las marcas pueden superponerse visualmente. La unión de la máscara puede contener varios runs separados; el gráfico los presenta posición por posición. Los archivos `region_verde.txt` y `region_verde.py` de esta copia solo informan `min(máscara)..max(máscara)`, por lo que pueden abarcar huecos y no enumeran los runs separados. Además, suman 1 al valor mínimo/máximo aunque `sirna_position` ya es 1-based y cortan la secuencia con `seq[ini:fin+1]`; no deben usarse como coordenadas canónicas. El explorador transforma la máscara en runs contiguos, conserva coordenadas 1-based inclusivas y mantiene aparte las alineaciones off-target.

## Datos por posición y unidades

Los resultados originales se generan en `siFi21--1.2.3/sifi_pipeline.py:SifiPipeline.data_to_json`. Los campos relevantes son:

| Campo JSON | Interpretación para la exploración |
|---|---|
| `sirna_position` | Inicio de ventana siRNA en query, 1-based |
| `sirna_sequence` | Secuencia de esa ventana; su longitud da el intervalo 1-based inclusivo del target |
| `is_off_target` | false: fila main-target; true: hit off-target; null: modo/contexto no informado |
| `is_efficient` | Etiqueta booleana original de si-Fi o null; no es probabilidad |
| `strand_selection` | Resultado booleano de criterio de selección de hebra, si fue calculado |
| `end_stability` | Resultado booleano del criterio de estabilidad de extremos, si se calculó |
| `target_site_accessibility` | Resultado booleano de umbral de accesibilidad |
| `accessibility_value` | Valor numérico extraído de la salida RNAplfold para la ventana/columna elegida |

En el JSON no se guardan energías numéricas, un score continuo de eficacia ni la máscara verde. `thermo_effcicient` tampoco se exporta. Las reglas terminales pueden alterar `is_efficient`; los componentes disponibles no permiten recalcularlo fielmente.

Las métricas originales se asignan a cada siRNA identificado por nombre, posición y secuencia. Un siRNA se incluye en el resumen de una ventana de comparación solo si **todo su intervalo cabe dentro de esa ventana**. Ventanas más cortas que un siRNA pueden no contener observaciones. Cada siRNA cuenta una vez aunque tenga varias filas de hit. Valores discrepantes del mismo siRNA no se promedian: quedan null/conflicto y se advierte. Se calculan conteos de true/false/null para etiquetas y filtros; `accessibility_value` recibe n, media, mediana, mínimo y máximo sobre observaciones numéricas disponibles. La proporción `is_efficient true / etiquetas true+false` es un resumen de etiquetas importadas, nunca una probabilidad ni una nueva predicción.

JSON si-Fi no conserva parámetros ni `sirna_size`. Se usa `sifi.run.parameters.sirna_size` si está en el manifiesto; de otro modo se infiere solo si todas las secuencias tienen una longitud uniforme, y se documenta esa inferencia. El usuario puede fijar `--sirna-size`. Si longitudes o posiciones son incompatibles con el candidato, la exploración se detiene. Métricas null permanecen ausentes; no se rellenan ni interpolan.

## Off-target asociado a cada ventana

El explorador examina las alineaciones descriptivas de `offtarget_evaluations[*].criterion.alignments`. Asocia a una ventana las filas cuyo intervalo sobre la query se superpone con ella, por organismo, y conserva estado de ejecución/criterio, ID de referencia, hebra, posiciones, longitud, mismatches, tipo, estado de coincidencia exacta, base/fingerprint y hash de resultado.

Cuando el importador dispone de inicio y longitud de query pero no fin explícito, calcula el fin inclusivo y marca `query_interval_derived_from_start_and_length=true`. No crea una alineación ni infiere que sea exacta. Una fila sin coordenadas query válidas no se atribuye a una ventana. Un criterio global de especie se conserva en cada contexto; un subconjunto local de ventanas no se usa para declarar un nuevo negativo o seguridad.

## Ejecutar, comparar y seleccionar

Primero debe existir un `report.json` que contenga resultados posicionales RNAi Designer. Los informes producidos solo por `analyze` sin `--sifi` dejan RNAi Designer como no evaluado y no tienen zonas verdes.

```powershell
python -m rnai_designer windows --report output/candidato/report.json --window-length 100 --step 10 --out output/candidato/windows-100
```

`window-length` y `step` están en nucleótidos. Se emiten únicamente ventanas completas contenidas en un run verde, empezando desde su extremo 1-based. `windows.csv` compara todos los resúmenes; `windows.svg` y `windows.png` grafican media y rango observado de accesibilidad y conteos de etiquetas. Una nueva corrida debe usar un directorio de salida vacío para no reemplazar resultados previos.

Para consultar la secuencia, siRNAs, métricas originales y filas off-target enlazadas a un inicio elegido de `windows.csv`:

```powershell
python -m rnai_designer windows --report output/candidato/report.json --window-length 100 --step 10 --select-start 301 --out output/candidato/window-301
```

La ventana seleccionada queda en `selected_window.json`; los solapamientos off-target por todas las ventanas están en `window_offtargets.csv`. `window_analysis.json` guarda parámetros, coordenadas, runs verdes, advertencias, hash del reporte fuente y el criterio off-target global. No modifica el reporte original ni las predicciones de si-Fi.
