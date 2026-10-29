# Optimización de rendimiento de rnai-designer 0.3.0

Se conserva la enumeración exhaustiva de 180 a 250 nt, paso de 1 nt, las regiones
verdes cortas, las coordenadas 1-based inclusivas, los puntajes y su orden de
desempate. No se modificaron si-Fi original, Conda legado, Bowtie, RNAplfold ni
las bases instaladas. Las modificaciones están dentro de rnai-designer.

## Medición y límites

`benchmarks/reference/` contiene una copia de los módulos anteriores a la
optimización; `benchmarks/results/environment.json` identifica sus hashes,
Python, sistema y procesador. `benchmarks/suite.py` ejecuta una sola versión y
un solo tamaño por proceso, secuencialmente. Los casos usan secuencias ACGT
generadas con semilla 817, una región verde completa y alineamientos
**explícitamente sintéticos**, uno cada diez posiciones. No se utilizan bases
biológicas ni se presentan los alineamientos fabricados como resultados reales.

Los tiempos de la tabla corresponden a validación del reporte, regiones verdes,
enumeración, métricas, puntuación, elección, CSV, JSON de selección y gráficos.
Excluyen arranque del intérprete e importaciones, preparación de los archivos del
benchmark y los ejecutables científicos externos. Son una observación por caso,
no una distribución estadística ni una garantía para cualquier equipo.

El pico de memoria es `PeakWorkingSetSize` de Windows, por proceso nuevo. Incluye
intérprete, bibliotecas, tablas y gráficos. No es simplemente el tamaño de las
estructuras Python. Las mediciones opcionales con `tracemalloc` y `cProfile` se
guardan por separado: su sobrecarga no se mezcla con los tiempos principales.

<!-- BENCHMARK_TABLE_START -->
| nt | Ventanas | Total antes → después (s) | Factor | Cálculo antes → después (s) | Exportación antes → después (s) | Pico RSS antes → después (MiB) | CSV antes → después (MB) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | 20,306 | 37.17 → 2.04 | 18.23× | 29.94 → 0.99 | 7.23 → 1.05 | 395.3 → 90.8 | 291.15 → 9.69 |
| 1000 | 55,806 | 117.80 → 4.61 | 25.58× | 90.68 → 2.94 | 27.11 → 1.67 | 944.0 → 118.3 | 806.44 → 26.66 |
| 2500 | 162,306 | 390.89 → 13.62 | 28.70× | 284.17 → 9.36 | 106.72 → 4.26 | 2251.4 → 173.2 | 2364.33 → 77.96 |
| 5000 | 339,806 | 1620.04 → 28.13 | 57.60× | 1029.51 → 19.02 | 590.54 → 9.10 | 2570.7 → 274.5 | 4962.40 → 163.58 |

En 5000 nt, el pico RSS disminuyó 89.32 %. Exportación incluye CSV, JSON de selección y figuras. MiB = 2²⁰ bytes; MB = 10⁶ bytes.

Durante la exportación de referencia de 5000 nt se observó además un máximo acumulado de 6.01 GB de memoria comprometida, registrado en `reference_5000_memory_observation.json`. Es un límite inferior del pico final, distinto del RSS. La ejecución optimizada final registró 0.54 GB de compromiso máximo.
<!-- BENCHMARK_TABLE_END -->

`benchmarks/results/comparison.json` conserva los tiempos sin redondear, tamaño
total de archivos, tamaño del CSV, memoria y número de filas verificadas. La
suite compara **cada fila** del CSV anterior con el nuevo, aplanando únicamente
sus componentes estructurados, y compara la ventana destacada completa salvo
la procedencia del archivo de entrada de cada benchmark. No se usan tolerancias
para ocultar diferencias numéricas.
La última mejora de almacenamiento se volvió a medir y se contrastó byte por
byte con esos CSV ya verificados, además de comparar nuevamente el ganador.
Las mediciones optimizadas anteriores se conservan en el JSON de comparación;
la tabla usa la implementación final, no mezcla tiempos con y sin instrumentación.

`benchmarks/stages.py` mide además preparación/validación, importación del JSON
si-Fi, resumen descriptivo, máscara verde, lectura de alineamientos Bowtie,
evaluación off-target, serialización JSON, reportes integrados con gráficos y
consolidación tabular de dos referencias. Los tiempos están en
`benchmarks/results/stages.json`. La consolidación aislada mide escritura
tabular; las pruebas del workflow verifican la integración real de múltiples
referencias y sus fallos parciales.

El entorno Conda legado no está disponible. Su prueba de integración se omite
con una razón explícita, y su tiempo figura como `null`, no como cero. Los
ejecutables **Bowtie 0.11.3 y RNAplfold 2.1.5** se encontraron posteriormente en
`siFi21--1.2.3/ToCopy`, fuera de las rutas configuradas por defecto.
`benchmarks/native_tools.py` los ejecuta directamente con FASTA e índices
temporales sintéticos. `native_tools.json` registra sus tiempos, construcción
de índices, búsqueda en tres especies y reutilización entre dos referencias,
con identidad de criterios comprobada. No se afirma una aceleración del motor
científico en sí ni un tiempo del pipeline completo a través de Conda.
Cada workflow genera
`performance.json`; la etapa headless genera `headless_performance.json`.

## Cuellos de botella encontrados

El perfil de referencia de 500 nt registra unos 272 millones de llamadas.
`_window_metrics` y sus agregados dominan el cálculo; `_gc_percent` se invoca
aproximadamente 3,95 millones de veces para solo 20.306 ventanas. Se recalculaban
clasificaciones, GC, cobertura, diversidad y estadísticas sobre los mismos
siRNAs solapados.
La puntuación se calcula dentro de `_window_metrics`; su tiempo está incluido
en ese tramo del perfil, junto con los agregados de cada ventana.

Cada ventana retenía copias de los alineamientos off-target solapados, incluyendo
procedencia. El exportador ejecutaba `json.dumps` sobre esas listas y varios
diccionarios por fila. El CSV de referencia de 2500 nt llegó a unos 2,36 GB.
El traceback observado durante `writer.writerows` era consistente con ese
trabajo repetido; no era el único costo.

También se construía un análisis completo de ventanas mínimas para validar la
entrada, antes de enumerar todas las longitudes. La consolidación de referencias
retenía nuevamente todas las ventanas. Los CSV integrados repetían las
advertencias y los parámetros completos de un dataset en cada registro: con
muchos datos faltantes, esto añadía otro crecimiento cuadrático de salida.

Los fingerprints leían íntegramente todos los índices en cada búsqueda. El mismo
gen se buscaba otra vez en las mismas tres especies por cada main target.

## Implementación

- **Sumas prefijas enteras** para etiquetas, completitud, conflictos, coincidencia
  con referencia, GC válido/objetivo, posiciones representadas y cobertura.
  Se validan una vez las longitudes uniformes de los siRNAs.
- **Ventana deslizante** para diversidad mediante `Counter` y para accesibilidad
  ordenada mediante `bisect`. Se reinicia al cambiar de rango o retroceder las
  coordenadas. No se interpolan datos ausentes.
- **Media de accesibilidad en el orden original**: se conserva `sum(access)`
  sobre la misma subsecuencia. No se resta una suma prefija flotante porque
  podría cambiar el redondeo y eventualmente un desempate. La fórmula de la
  mediana y los puntajes mantiene el orden aritmético anterior.
- **Índices de intervalos off-target**: para contar solapamientos se usa
  `cantidad(inicio <= fin_ventana) - cantidad(fin < inicio_ventana)` sobre arrays
  ordenados, por especie y categoría. El detalle se materializa solo cuando se
  solicita o para la ventana destacada.
- **Máscara verde con bitsets enteros** que conservan las instantáneas y el orden
  del helper original. No se reemplaza la semántica histórica por una resta
  global distinta. Pruebas con filas mezcladas, varios hits y `null` contrastan
  el resultado con la copia anterior.
- **Tabla columnar sin dependencias nuevas**: `CompactWindows` usa `array` de la
  biblioteca estándar. No conserva un diccionario ni listas de alineamientos
  por ventana. Los textos repetidos usan categorías; el identificador se deriva
  de región y coordenadas. Se mantiene acceso por índice e iteración como
  diccionarios, que se reconstruyen a demanda. No debe convertirse toda la
  tabla a `list` para exportarla.
- **Ganador y cantidad de empates en una pasada**, con la misma tupla de
  comparación. No se construye la lista adicional de todas las filas rankeables.
- **Esquema de columnas preparado una sola vez**: se anexan valores directamente
  a los arrays. No se reconstruyen diccionarios aplanados ni categorías vacías
  en cada fila. Las columnas numéricas mixtas promueven a flotante para impedir
  truncamientos de decimales; sus valores originales siguen en el reporte fuente.
- **CSV incremental de escalares**: el exportador recorre columnas mediante
  iteradores, sin `json.dumps` por ventana. Los CSV de múltiples referencias
  usan un archivo temporal de disco y liberan la selección de cada referencia
  antes de seguir con la siguiente.
- **JSON incremental** mediante `json.dump` para el reporte integrado, evitando
  crear una cadena adicional del tamaño del reporte en memoria.
- **Gráficos con segmentos agrupados** mediante `LineCollection`: se evita un
  objeto gráfico por alineamiento. Se conserva el capstyle y se comprobó
  identidad de los píxeles PNG con la implementación anterior.

SQLite y formatos columnares externos no fueron necesarios: las columnas
compactas permiten acceso aleatorio compatible con el explorador y requieren
solo la biblioteca estándar. SQLite añadiría índices, transacciones y costo de
inserción para cientos de miles de filas. Se puede evaluar como opción de
almacenamiento en disco para secuencias mucho mayores; no se afirma un benchmark
de SQLite que no se haya ejecutado.

## Complejidad

`N`: siRNAs únicos; `R`: filas originales; `W`: ventanas; `K`: siRNAs por
ventana; `H`: alineamientos; `S`: nucleótidos del gen; `B`: bytes de índices.
El rango predeterminado mantiene `K <= 230` si hay una observación por inicio;
puede haber más observaciones distintas en una misma posición.

| Operación | Antes | Después |
|---|---|---|
| Flags, calidad y GC por ventana | O(K), repetidas pasadas | O(1) por agregado, tras O(N + S) |
| GC de cada siRNA | O(W × K × 21) | O(N × 21), con caché acotada de función pura |
| Diversidad | O(K) por ventana | Actualizaciones al avanzar, normalmente O(1) |
| Mediana | O(K log K) por ventana | O(K) por actualización de lista ordenada; lectura O(1) |
| Media de accesibilidad | O(K) | O(K), deliberadamente por igualdad aritmética |
| Cobertura de intervalos uniformes | O(K log K) | O(1), tras O(N) |
| Conteos off-target | O(W × H) | O(H log H + W log H), categorías constantes |
| Memoria de ventanas | O(W × campos + copias de hits solapados) | O(W × campos escalares) |
| Exportación de ventanas | O(W + bytes de detalles repetidos) | O(W × campos escalares), incremental |
| Consolidación en memoria | O(suma de ventanas de referencias) | O(una referencia), más spool de disco |
| Hashes repetidos de índices | O(B) cada vez | O(B) inicial; O(número de archivos) sin cambios |

La enumeración sigue siendo O(W). Se eliminó el conjunto auxiliar de intervalos
vistos porque las regiones reconstruidas son disjuntas y cada longitud/inicio
se genera exactamente una vez.

## Caché y reutilización

1. **Hashes por proceso**: SHA-256 se calcula en bloques de 1 MiB. Se reutiliza
   únicamente si ruta resuelta, dispositivo, identidad del archivo, tamaño,
   mtime y ctime no cambiaron. Se verifica también que el archivo no cambie
   durante la lectura. No se guarda una confianza persistente basada solo en
   fechas entre ejecuciones. Las bases deben permanecer estables durante un
   análisis; una modificación externa que falsee todos los metadatos queda
   fuera de ese contrato. Reiniciar el proceso fuerza la lectura completa.
2. **Búsquedas off-target entre referencias**: una sesión limitada al workflow
   multirreferencia reutiliza únicamente búsquedas terminadas. La clave incluye
   candidato completo y marca de simulación, ejecutable/dependencias, hashes y
   procedencia de bases y parámetros. La evidencia se verifica y copia a la
   nueva referencia; `reuse_manifest.json` identifica la ejecución original.
   Se comparten en memoria los registros científicos, sin copiarlos por ventana.
   Los errores no se guardan como resultados completados. Una nueva sesión
   comienza sin este caché. Se retienen como máximo cuatro búsquedas mediante
   LRU para no acumular todas las ventanas destacadas de todas las referencias;
   la búsqueda compartida del gen se mantiene reciente al reutilizarse.
3. **Caché headless opcional**: `workflow --cache-dir <directorio>` permite
   reutilizar los JSON originales, la máscara exportada y el archivo crudo de
   alineamientos `.bowtie.tmp`, sin descartar evidencia. Incluye candidato
   completo, simulación, hash del FASTA, main target, parámetros, versión de
   rnai-designer, información del runtime, intérprete, metadatos de paquetes
   Conda, código original de si-Fi, adaptador, herramientas/dependencias e
   índices. Las ubicaciones de salida no forman parte de la identidad del
   FASTA. Los artefactos se verifican por contenido antes de restaurarlos.
   Una entrada sin `complete.json`, corrupta o con clave diferente es un miss.
   Las entradas se publican atómicamente. No se escribe dentro de las bases.

No se utiliza el resultado del benchmark como entrada real. Las pruebas de
caché separan explícitamente candidatos simulados y reales; la ventana derivada
conserva ahora la marca de simulación de su gen.

La caché headless evita repetir el motor completo cuando coinciden sus entradas;
la enumeración inicial de targets aún puede invocar Bowtie. No se asegura una
reanudación de cada operación del pipeline: se pueden reutilizar artefactos
verificados desde una salida nueva, y los main targets se siguen procesando
independientemente tras fallos de otras referencias.

## Escrituras, progreso y migración

Los reportes integrados, exploración y selección se escriben en directorios
hermanos `.partial` y se renombran al finalizar todos sus archivos y figuras.
Si el reporte integrado comparte una carpeta existente con evidencia, se
reemplazan únicamente sus archivos generados, conservando los demás. Se publica
`report.json` al final y se revierten los reemplazos si falla la publicación.
El headless escribe sus JSON en un temporal y verifica su estructura antes de
publicarlos. Los resúmenes y la consolidación usan archivos temporales con
reemplazo atómico. Un fallo de escritura o `KeyboardInterrupt` no publica la
selección ni reemplaza un resumen válido por un archivo truncado. Un cierre
forzado puede dejar temporales `.partial`, que no son resultados completados.
Esto no equivale a una garantía de durabilidad frente a pérdida de energía.
Los directorios publicados heredan los permisos de su carpeta de destino;
se evita publicar el ACL restrictivo de `mkdtemp` en Windows recientes.

Se informa progreso por etapas, y se registran tiempos reales en
`performance.json` sin alterar el resumen científico. No se emite una actualización
por cada ventana, para evitar otro costo de E/S.

**CSV 2.0** (el JSON integrado conserva su esquema y se añade
`csv_schema_version`):

- `all_windows.csv` y `all_references_windows.csv` mantienen las columnas
  escalares científicas. Los componentes usan nombres
  `quality_components.<nombre>` y `efficacy_priority_components.<nombre>`.
  Se elimina `offtarget_matches`; `offtarget_by_species` ya está representado en
  las columnas escalares por especie. El detalle permanece en el reporte de
  entrada y en la ventana destacada.
- `records.csv` conserva los campos científicos y `record_json` para los
  consumidores existentes. Sustituye `provenance_json`, `parameters_json` y
  `warnings_json` repetidos por `dataset_ref` (JSON Pointer dentro de
  `report.json`) y `record_index` (índice desde cero). Los registros originales
  tienen un repositorio autoritativo en el JSON integrado; se conserva su
  proyección CSV y la copia archivada de entrada para compatibilidad/trazabilidad.
- Las columnas de métricas que mezclan enteros y decimales se representan como
  flotantes sin truncar sus valores. Los tipos originales permanecen en el JSON
  científico de entrada.
- `evaluations.csv` sustituye la copia completa `evaluation_json` por
  `evaluation_ref`, también un JSON Pointer. El detalle no se pierde.
- La API Python devuelve una secuencia `CompactWindows`, no una lista mutable.
  Se conserva `len`, indexación, slices e iteración de diccionarios. Los slices
  materializan sus filas. La tabla completa no se pasa directamente a
  `json.dumps`; para exportaciones se usan los escritores incrementales.

Recuperar detalle de una ventana concreta:

```python
from rnai_designer.data_models import read_json
from rnai_designer.window_details import resolve_window_details

report = read_json('integrated/report.json')
detail = resolve_window_details(report, 'green-001:100-279')
# detail contiene secuencia, siRNAs con sus filas originales y off-targets.
```

El resultado destacado sigue en `selected_window.json`, con su filtro y su
política original de liberación de secuencia. El filtro no se aplica para
cambiar el ganador ni se selecciona otro silenciosamente si falla.

## Verificación y reproducción

```powershell
.venv/Scripts/python.exe -m pytest -q -rs
.venv/Scripts/python.exe benchmarks/suite.py
.venv/Scripts/python.exe benchmarks/stages.py
.venv/Scripts/python.exe benchmarks/native_tools.py
.venv/Scripts/python.exe benchmarks/benchmark.py --mode before --lengths 500 --profile
.venv/Scripts/python.exe benchmarks/benchmark.py --mode after --lengths 500 --profile
.venv/Scripts/python.exe benchmarks/benchmark.py --mode before --lengths 500 --trace
.venv/Scripts/python.exe benchmarks/benchmark.py --mode after --lengths 500 --trace
.venv/Scripts/python.exe benchmarks/parallel_probe.py
```

Las mediciones de referencia pueden producir varios GB de CSV; no deben
ejecutarse sobre bases biológicas. Los directorios de exportación del benchmark
se ignoran en Git; se conservan los JSON de medición y perfiles de texto.

Las pruebas diferenciales usan igualdad numérica exacta, sin tolerancia flotante:
flags, cobertura, GC, media/mediana, puntajes, desempates, ventana destacada,
conteos por especie, duplicados, conflictos, datos ausentes, bases ambiguas,
regiones cortas e incluso regiones menores al tamaño de siRNA. Se conservan las
pruebas existentes y se añaden fallos de escritura, interrupciones, caché
obsoleta/incompleta/corrupta, separación de simulación y múltiples referencias.

<!-- VALIDATION_START -->
120 passed, 1 skipped in 54.27s

La omisión corresponde exclusivamente al headless Conda legado. La prueba de Bowtie real usa los ejecutables empaquetados y un índice sintético temporal.

- `tracemalloc`, before, 500 nt: pico Python 310.16 MiB; tiempo instrumentado 290.51 s (excluido de la tabla principal).
- `tracemalloc`, after, 500 nt: pico Python 13.42 MiB; tiempo instrumentado 32.69 s (excluido de la tabla principal).

Sonda Windows, dos referencias sintéticas de 500 nt: secuencial 6.23 s; spawn, cálculo y transferencia 4.24 s; factor 1.47×. Se verificaron todas las filas y el ganador. No se habilitó paralelismo del producto: esta sonda no valida el rendimiento ni el aislamiento de Conda/RNAplfold simultáneos.

| nt | Bowtie tres especies, primera búsqueda (s) | Repetición anterior (s) | Reutilización verificada (s) | RNAplfold directo (s) |
|---:|---:|---:|---:|---:|
| 500 | 0.926 | 0.660 | 0.326 | 0.198 |
| 1000 | 1.513 | 0.689 | 0.375 | 0.637 |
| 2500 | 1.012 | 0.774 | 0.228 | 1.318 |
| 5000 | 1.124 | 0.954 | 0.340 | 2.438 |

La primera búsqueda incluye fingerprints, versión, ejecución y parseo. El JSON conserva también el tiempo aislado de `execute` y el número de invocaciones: tres al buscar y cero al reutilizar. La caché puede añadir costo a la primera búsqueda; su beneficio está en las repeticiones. Las herramientas nativas no se aceleraron ni se modificaron: se evitó repetir búsquedas completas. Se verificó igualdad de los criterios científicos en ambos modos. Los índices son pequeños y sintéticos; estos tiempos no estiman búsquedas en las bases biológicas instaladas.
<!-- VALIDATION_END -->

## Compromisos y siguientes mejoras

El volumen tabular sigue siendo proporcional a todas las ventanas requeridas;
no se reduce ni el rango ni la exhaustividad. La memoria no es O(1): la API y
el acceso aleatorio conservan columnas para todas las ventanas. El detalle
extenso sí se obtiene a demanda. El spool de referencias intercambia memoria
por E/S secuencial de disco.

Se mantiene ejecución científica secuencial. Las herramientas legacy usan
directorios de trabajo compartidos y no se probó su concurrencia en esta
máquina. La sonda de multiprocessing mide únicamente el cálculo Python y el
traslado de tablas, incluyendo spawn de Windows; no habilita paralelismo en el
producto. Un beneficio aislado no demuestra un beneficio del pipeline completo.

Próximos pasos posibles: medir el headless completo con Conda instalado,
repetir benchmarks para intervalos de confianza, almacenamiento
columnar en disco para genes mucho mayores, mediana con estructura de orden
logarítmica si se amplía mucho el rango, y reanudación granular de búsquedas
externas entre invocaciones. Cambiar la media por restas de sumas prefijas
flotantes requeriría demostrar que no modifica puntajes ni desempates; no se
incluyó esa aproximación.
