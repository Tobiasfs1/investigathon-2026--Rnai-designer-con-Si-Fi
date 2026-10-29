# rnai-designer

La versión 0.3.0 incorpora tablas compactas, exportación incremental y caché
headless opcional (`workflow --cache-dir <directorio>`). Los cambios de CSV 2.0,
benchmarks, conservación de resultados y tiempos por etapa están documentados en
[PERFORMANCE_OPTIMIZATION.md](docs/PERFORMANCE_OPTIMIZATION.md).

Módulo Python 3.10+ para integrar RNAi Designer con su GUI e importar sus métricas, y automatizar búsquedas de secuencias off-target con el Bowtie instalado por si-Fi. Conserva el puntaje upstream. Las búsquedas por Bowtie no calculan eficiencia ni accesibilidad. Los criterios computacionales de coincidencia son condicionales a las bases consultadas y no prueban seguridad biológica ni ambiental.

## Ejecutar el ejemplo

Desde `rnai-designer/`, instalar la dependencia de gráficos y ejecutar:

```powershell
python -m pip install -r requirements.txt
python -m rnai_designer --candidate examples/candidate.json --sifi examples/sifi_result.json --run examples/run.json --offtargets examples/offtargets.json --out output/demo
```

Abrir `output/demo/efficiency.svg` en el navegador. Revisar `report.json`, `records.csv`, `positions.csv` y `evaluations.csv`. **Todos los datos de examples son FICTICIOS**, escritos manualmente para probar el formato: no representan genes, experimentos ni ejecuciones reales de si-Fi. Incluyen valores faltantes intencionales. Los ejemplos de resultados contra humano, maíz y abeja tampoco son evaluaciones reales.

Configuración automática en Windows: `%LOCALAPPDATA%\siFi2015\Databases`; se puede cambiar con `--config`, `--db-dir` o `SIFI_DB_DIR`. Verificar herramientas e índices:

```powershell
python -m rnai_designer doctor --config examples/config.windows.json
python -m rnai_designer databases --config examples/config.windows.json --inspect
```

Para ejecutar las búsquedas off-target automáticas contra las tres bases configuradas usando un candidato real y una salida nueva:

```powershell
python -m rnai_designer analyze --candidate ruta/candidato.json --config examples/config.windows.json --out output/candidato-2026-10-27
```

El análisis envía cada ventana 21 nt a Bowtie con cero mismatches, sin gaps y ambas hebras. Sus reportes no incluyen las métricas de eficiencia de RNAi Designer. No ejecutar fixtures sintéticos contra bases biológicas: usar el ID y la secuencia autorizados por el equipo. Un directorio de salida existente con archivos se rechaza para preservar resultados anteriores.

Para usar el modo manual, preparar un handoff y abrir su FASTA en la GUI. Si el entorno legado esta disponible, se puede ejecutar el motor original en headless mediante el comando workflow documentado mas abajo:

```powershell
python -m rnai_designer prepare --candidate ruta/candidato.json --out output/candidato-designer
```

Después de ejecutar el diseñador y completar `output/candidato-designer/run.json`, recuperar un JSON temporal compatible:

```powershell
python -m rnai_designer collect --candidate ruta/candidato.json --run output/candidato-designer/run.json --search-dir $env:TEMP --out output/candidato-resultados
```

Si aparecen varios archivos JSON compatibles, el importador se detiene para que se indique `--sifi` explícitamente. Detalles y límites en [docs/SIFI_INTEGRATION.md](docs/SIFI_INTEGRATION.md).

## Explorar ventanas verdes

Cuando el reporte ya contenga resultados de RNAi Designer, `windows` reconstruye las regiones main-target y compara ventanas enteras de longitud configurable:

```powershell
python -m rnai_designer windows --report output/candidato-resultados/report.json --window-length 100 --step 10 --out output/candidato/windows-100
```

Seleccionar un inicio 1-based que aparezca en `windows.csv` y repetir con una carpeta nueva y `--select-start 301` para exportar el detalle. Las tablas/figuras comparan métricas descriptivas originales, nunca nuevas predicciones. Definición del verde, coordenadas y formatos: [docs/WINDOW_EXPLORATION.md](docs/WINDOW_EXPLORATION.md).

Opcional: `python -m pip install -e .` habilita el comando `rnai-designer`. Para pruebas:

```powershell
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

`requirements-tested.txt` fija las versiones usadas para validar esta entrega en Windows con Python 3.14.10. Para reproducir ese entorno: `python -m pip install -r requirements-tested.txt`. En otras versiones de Python usar los rangos de requirements-dev.txt; no se ejecutó una matriz de compatibilidad. Los reportes conservan versión del módulo y hashes de las entradas para rastrear una ejecución.

## Importar resultados reales

1. Ejecutar si-Fi por separado en su entorno legado. Usar una query por ejecución y registrar parámetros y base.
2. Copiar el JSON temporal generado por `SifiPipeline.run_pipeline` a una carpeta estable. No sirve el CSV resumen. La GUI inspeccionada no ofrece exportación explícita del JSON: se necesita identificarlo en TEMP y verificar query_name y contenido. No elegirlo solo por fecha.
3. Crear candidate.json con referencia exacta, procedencia, versión, ID del homólogo y prioridad upstream. Crear run.json según examples, reemplazando todos los datos ficticios. Dejar desconocidos explícitos; no inventar parámetros.
4. Ejecutar la CLI con esas rutas. Sin `--offtargets`, las tres especies figuran no evaluadas. Para incorporarlas, obtener resultados externos contra bases documentadas por especie y crear el manifiesto de [docs/FORMATS.md](docs/FORMATS.md).

Los errores de lectura o formato devuelven código 2; las advertencias quedan en JSON. Los resultados de esta copia de si-Fi omiten hits de hebra negativa y tienen una rama sin hits defectuosa: no bastan para declarar seguridad. El SVG muestra accesibilidad y etiquetas del software, no un score continuo de eficacia.

## Arquitectura y límites

Módulos independientes: `data_models`, `sifi_adapter`, `efficiency_analysis`, `offtarget_reports`, `visualization` y `reporting`; la CLI los coordina. Matplotlib genera SVG y PNG sin interfaz gráfica. Contratos: [docs/FORMATS.md](docs/FORMATS.md). Auditoría con archivos y funciones originales: [docs/SIFI_ANALYSIS.md](docs/SIFI_ANALYSIS.md).

Pendiente validar un candidato real con procedencia y bases versionadas para humano, maiz y abeja. El motor original ya se pudo ejecutar en headless sobre una referencia sintetica temporal; esto no es validacion biologica. No se implementa ranking de eficacia. El codigo original y las bases permanecen intactos.

### Selección automática de ventanas

Evalua todas las ventanas 180-250 nt con paso de 1 nt y agrega una ventana exacta para cada run verde mas corto. Prioriza mediante la heuristica exploratoria de GC por siRNA, diversidad y etiquetas originales `is_efficient`; conserva la calidad de evidencia por separado.

~~~powershell
python -m rnai_designer select-windows --report output/candidato-resultados/report.json --out output/candidato/selection
~~~

Genera selected_window.json, all_windows.csv, selection_report.json y una visualizacion del gen. La secuencia solo se publica si el filtro Bowtie de la ventana registra cero coincidencias exactas de 21 nt en humano, maiz y abeja con las tres corridas completas. La calidad de evidencia queda como auditoria independiente.

## Flujo integrado

Con `--run-sifi-headless`, workflow detecta todos los main targets hallados por Bowtie y procesa cada referencia por separado. Cada referencia conserva su JSON, mascara verde, coordenadas, off-targets y seleccion en `references/`; el reporte consolidado y la tabla conjunta quedan en la raiz.

~~~powershell
python -m rnai_designer workflow --gene-id DM_gen --sequence-file gen.fasta --reference-source repositorio --reference-version accesion --sifi ruta/sifi.json --run ruta/run.json --config examples/config.windows.json --out output/DM_gen
~~~

Para iniciar primero con handoff, usar `workflow` sin `--sifi`. Para que el mismo proceso espere mientras se opera la GUI legada, agregar `--wait-for-sifi`; el programa detectará el JSON temporal y continuará. La GUI de si-Fi sigue requiriendo la selección manual de main targets. Detalles y límites de ejecución: [docs/SIFI_AUTOMATION.md](docs/SIFI_AUTOMATION.md).


## Ejecuci?n integrada con el si-Fi legado

En Windows, `doctor` detecta Conda y el entorno antiguo si esta disponible. Para ejecutar todos los main targets con una sola entrada, usar `workflow --run-sifi-headless` sin `--main-target`. Para limitarse a una referencia, indicar `--main-target ID`; `--list-targets-only` enumera las referencias disponibles.

```powershell
python -m rnai_designer workflow --gene-id ATPase --sequence-file inputs/ATPase.fasta --out output/ATPase_multi --run-sifi-headless
```

El adaptador requiere el entorno Python 2/PyQt4 existente, ademas de Bowtie y RNAplfold. No instala dependencias ni modifica si-Fi o sus bases. Guarda el JSON original, la mascara verde y la procedencia. Workflow publica la secuencia priorizada solo si las tres busquedas Bowtie de la ventana concluyen completas y sin hits exactos de 21 nt; la prioridad es heuristica y no eficacia validada. Consultar docs/SIFI_AUTOMATION.md.


## Portabilidad entre computadoras Windows

Ver [instalacion portable, configuracion y doctor](docs/portabilidad.md).


El comando `workflow` muestra por defecto un resumen tabular y la secuencia de la
ventana priorizada, junto con el estado computacional off-target. La secuencia
solo se muestra cuando el filtro existente autoriza su presentacion.
El resumen compacto es el predeterminado (unas 23 lineas mas la secuencia,
en bloques de 60 bases). Usar `--verbose-report` para el detalle por referencia.
El formato de consola no altera los reportes JSON y CSV.
`--length-penalty-weight 30` configura el peso (finito y no negativo) de
`length_penalty = peso * max(0, 1 - length_nt / 180)`.
La clasificacion usa `adjusted_priority_score = max(0, efficacy_priority_score - length_penalty)`;
el puntaje original y todas las metricas cientificas se conservan.
Sin puntaje original, el ajustado tambien es nulo y la ventana no es elegible.
Las clases `short_green_region` y `conventional_180_250nt` siguen separadas;
el primer puesto global es exploratorio y no establece equivalencia biologica.
Los empates usan los componentes originales, calidad, coordenadas y referencia.
Los reportes antiguos se pueden mostrar con `console_output.print_workflow_summary`;
los campos nuevos ausentes se indican como no disponibles, sin inventar evidencia.


### Ejecucion simple

Desde la carpeta de `rnai-designer`:

```powershell
python -m rnai_designer MiGen "ATGCGTACGATCG..."
```

Reemplazar los puntos por la secuencia completa IUPAC. Se aceptan minusculas,
espacios y saltos de linea dentro del argumento entre comillas; se normalizan
con el mismo contrato de entrada del workflow avanzado.
Se requieren exactamente dos argumentos. No es necesario preparar archivos.

Se crean `inputs/MiGen_1.fasta` y `output/MiGen_1/` automaticamente, relativos
a la carpeta actual. Si el FASTA o la carpeta ya existen, se busca el siguiente
numero libre. El identificador base se sanea para Windows y se limita a 12
caracteres antes de agregar el correlativo; no se usan fechas ni UUID para estos
nombres. El nombre completo se conserva en `gene_id` y en los metadatos de
`candidate.json`. Se comprueba el limite de ruta de Windows reservando margen
para archivos internos; si la carpeta actual es demasiado larga, se informa
el error antes de crear los archivos.

La ejecucion delega al mismo `workflow --run-sifi-headless` con sus valores
predeterminados y configuracion existente, detecta los targets y conserva el
resumen compacto. Requiere las mismas herramientas y bases configuradas que
el workflow avanzado; no instala ni cambia el entorno. Si falla el analisis,
los archivos reservados permanecen para inspeccion y la siguiente ejecucion
usa otro numero. Los comandos `workflow`, `doctor` y
`workflow --resume output/MiGen_1 --run-sifi-headless` siguen disponibles.
Para opciones avanzadas, incluida `--verbose-report`, usar `workflow`.


### Analisis batch

```powershell
python -m rnai_designer batch genes.fasta
```

El FASTA contiene una entrada por gen. Se utiliza el encabezado completo como
nombre original. Se valida todo el archivo antes del analisis; las secuencias
invalidas se registran individualmente y no bloquean los otros genes. Una
estructura FASTA incorrecta (por ejemplo, secuencia sin encabezado) impide iniciar.
Los encabezados duplicados se conservan como entradas independientes, identificadas
por su posicion `input_index` en el archivo.

Se crea una carpeta corta nueva `output/batch_1/` (o el siguiente numero libre).
Dentro quedan los FASTA en `inputs/g1_1.fasta`, `inputs/g2_1.fasta`, etc., los
reportes completos en `output/g1_1/`, `output/g2_1/`, y un log por intento.
Cada gen ejecuta el mismo workflow headless y configuracion predeterminada que
el comando individual. La consola muestra progreso y un ranking final compacto.

La carpeta batch contiene:

- `resultados_genes.csv`: una fila por entrada, conservando orden y estado.
- `ranking_genes.csv`: primero las entradas evaluables, ordenadas por puntaje
  ajustado descendente; los empates conservan el orden de entrada. Las entradas
  sin puntaje ajustado quedan al final como `no_evaluable_result`, sin rango ni
  puntaje inventado. Un cero calculado sigue siendo un resultado evaluable.
- `resumen_batch.json`: entradas normalizadas, estados, errores, intentos,
  ubicaciones de los analisis, resultados y huellas de integridad.

Los CSV conservan la clase de ventana, suficiencia de datos, estado de ejecucion,
liberacion de secuencia y estados off-target por especie. Los resultados parciales
pueden tener ventanas evaluables y se identifican por su estado. El ranking es
heuristico y no constituye una prediccion validada de eficacia biologica.

Para reanudar, sin necesitar el FASTA original:

```powershell
python -m rnai_designer batch --resume output/batch_1
```

Tambien se admite `batch genes.fasta --resume output/batch_1`; en ese caso se
comprueba que el FASTA coincide con el original. Los genes completados con reporte
integro no se repiten. Los intentos fallidos o interrumpidos se reejecutan en otra
carpeta corta, conservando la evidencia anterior. Si un reporte final ya se publico
antes de una interrupcion del batch, se recupera sin repetir el analisis. Los genes
con entrada invalida permanecen registrados; para corregirlos se inicia otro batch.
Una ejecucion simultanea del mismo batch se rechaza mediante bloqueo del sistema.
Los reportes generales se actualizan atomicamente tras cada gen. Codigo de salida:
0 si todos completan, 2 si hay errores y 130 si se interrumpe con Ctrl+C.


### Seleccion final condicionada al filtro off-target

El workflow (incluidos el comando simple y batch) recorre las ventanas elegibles
en el mismo orden de puntaje ajustado y desempate existente. Si la primera no
cumple el criterio computacional, evalua la siguiente, hasta encontrar una
aprobada o agotar las candidatas. Los puntajes originales, penalizaciones,
metricas y clases de ventana no cambian.

La ganadora debe acreditar cero coincidencias exactas de 21 nt en humano, maiz y
abeja mediante el filtro existente. Un error tecnico, un resultado indeterminado
o una evaluacion no realizada no autorizan la seleccion: se registran y se
continua. No hay un limite artificial de intentos, por lo que el tiempo puede
crecer considerablemente si muchas ventanas fallan. Si ninguna acredita el
criterio, `selected_window` es null y no se publica secuencia ganadora.

`offtarget_selection` en los JSON de seleccion y pipeline registra cada intento,
su posicion original en el ranking, puntajes, filtro, motivos y carpeta de
resultados. La ganadora incluye `selected_priority_rank`. Las corridas de la
primera candidata se conservan en `offtarget_selected_window/`; las posteriores
en `offtarget_w2/`, `offtarget_w3/`, etc. El nombre de la primera carpeta se
mantiene por compatibilidad y no significa que haya sido aprobada.
La comparacion entre referencias elige entre sus ganadoras aprobadas. El ranking
batch excluye ventanas cuyo filtro no acredite aprobacion, tambien al leer
reportes anteriores; los reportes historicos no se reescriben. El criterio es
computacional y acotado, no una garantia general de seguridad biologica.


### Reutilizacion por gen (sin paralelismo)

Durante un workflow headless, incluida cada entrada batch, se utiliza un espacio
compartido temporal y privado. La enumeracion de targets y los GITV posteriores
reutilizan los alineamientos Bowtie de Dalbulus. RNAplfold y los registros base de
si-Fi se calculan una vez cuando se necesitan; para cada GITV se actualiza solo
`is_off_target`, la unica etiqueta de `data_to_json` dependiente del target. La
mascara verde sigue pasando por `general_helpers.get_target_data` original,
conservando el orden de filas y sus snapshots. Los resultados de eficiencia con
identicos argumentos se reutilizan dentro de la configuracion fija del pipeline.
No se cambia el codigo cientifico original ni se filtran GITV por numero de hits.

El espacio temporal se elimina al terminar o salir por excepcion; no se comparte
entre genes ni ejecuciones. Su binding comprueba secuencia, ID, parametros,
rutas y firmas de dependencias antes de reutilizar resultados; si cambian durante
el workflow, se rechaza la reutilizacion. Esto es independiente de `--cache-dir`.
Las evidencias de alineamientos y los JSON/mask individuales siguen exportandose.
Siguen existiendo el arranque Conda por referencia y las escrituras individuales;
no se promete un factor de aceleracion medido.

La seleccion segura prueba la mejor ventana conocida antes de ordenar las demas;
si hace falta continuar, usa columnas escalares de la tabla compacta y conserva
exactamente los desempates. La preparacion de siRNAs se comparte entre validacion
y seleccion; sin regiones verdes se evitan indices y calculos exclusivos de
ventanas, manteniendo los reportes, advertencias y evaluaciones del gen completo.
La cache off-target retiene la corrida del gen completo durante las referencias
sin que la expulsen las ventanas rechazadas; mantiene como maximo cuatro corridas
de ventanas adicionales. Tambien se evita volver a leer el JSON si-Fi que ya se
cargo al validar la exportacion, y se copian evidencias grandes sin cargarlas
enteras en memoria.
