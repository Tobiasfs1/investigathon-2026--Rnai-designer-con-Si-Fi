# Batch, ranking y reportes

[Índice](README.md) · [Criterios científicos](criterios-cientificos.md) · [Almacenamiento](rendimiento-y-limitaciones.md)

## Entrada y ejecución

Desde `rnai-designer/`, con Python 3 configurado:

```powershell
python -m rnai_designer batch genes.fasta
```

Cada encabezado completo del FASTA es el nombre original de una entrada. Se normalizan mayúsculas y espacios de la secuencia y se valida IUPAC. Encabezados duplicados se conservan como entradas independientes mediante `input_index`; no se fusionan por nombre.

Se valida el archivo antes de comenzar. Una entrada con secuencia inválida se registra como `invalid_input` y no bloquea las válidas. Una estructura ilegible, por ejemplo nucleótidos antes del primer encabezado, impide iniciar el batch.

No se requiere otra opción para elegir parámetros: cada gen utiliza el workflow headless y la configuración predeterminada vigente. El comando simple y el batch admiten menos opciones que `workflow`; no asumir que aceptan todos sus flags.

## Timeout por gen

Cada intento se ejecuta en un trabajador Python independiente, con límite total de **360 segundos (6 minutos)** por defecto. El contador comienza antes de preparar la salida del intento e incluye el lanzamiento/importación del trabajador, si-Fi21, las ventanas y las evaluaciones off-target. No es un límite separado por GITV ni sustituye los timeouts internos de las herramientas.

```powershell
python -m rnai_designer batch genes.fasta
python -m rnai_designer batch genes.fasta --gene-timeout 600
python -m rnai_designer batch genes.fasta --gene-timeout 0
python -m rnai_designer batch --resume output/batch_1 --gene-timeout 600
```

`--gene-timeout` acepta segundos finitos no negativos. **0 desactiva el límite**, manteniendo el trabajador independiente y la limpieza de procesos. En reanudación se conserva `gene_timeout_seconds` del resumen si no se pasa otra opción; para un batch anterior sin ese campo se aplica 360.

Al vencer el límite se cancela el árbol de procesos antes de continuar. En Windows se usa un Job Object, asignado antes de autorizar el trabajo; incluye sus descendientes y termina al cerrar el último handle. No se usa un hilo que deje el cálculo original activo. Fuera de Windows se utiliza un grupo de procesos. La cuenta final `elapsed_s` incluye el tiempo de cancelación/cierre y puede superar levemente el límite.

El estado del gen y del intento es **`incomplete_timeout`**, nunca aprobación ni rechazo científico. Los CSV mantienen sus columnas: el estado y error indican el timeout, el puntaje queda vacío y la entrada no participa como resultado evaluable en el ranking. `resumen_batch.json` conserva `elapsed_s`, `timeout_seconds`, `timeout_stage`, `timeout_reference_id` y, por intento, PID y confirmación de limpieza. La etapa y referencia son las últimas publicadas por los mensajes de progreso existentes; si el motor no publicó esa información, pueden no estar disponibles.

Los genes ya completados permanecen intactos. Un gen vencido se reintenta en otra carpeta con `--resume`; no se recupera como definitivo un reporte parcial aunque el archivo aislado diga `complete`. Los intentos nuevos requieren también constancia de finalización satisfactoria del trabajador dentro del límite. Se mantienen las reglas de recuperación anteriores para intentos históricos sin metadatos de trabajador.

Se agregan archivos pequeños junto al log de cada intento: `g1_1.worker.json` (supervisión), `g1_1.done.json` (constancia del trabajador si llegó a terminar), `g1_1.progress.json` (último contexto) y `g1_1.request.json` (entrada al trabajador). Un área temporal privada se elimina después de terminar sus procesos; se conservan las evidencias parciales en la carpeta del intento.

El soporte Windows requiere poder crear/asignar Job Objects; los jobs anidados son compatibles desde Windows 8. Si no se puede establecer la contención, no se inicia el análisis científico de ese trabajador. Si no se puede confirmar una cancelación, el batch se detiene con error en lugar de iniciar otro gen con procesos posiblemente activos. Referencia técnica: [Microsoft, Job Objects](https://learn.microsoft.com/en-us/windows/win32/procthread/job-objects).

El límite solo afecta al comando `batch`; no cambia el funcionamiento individual de `workflow`, las fórmulas, regiones o criterios off-target. Una ejecución que ya estaba iniciada no adquiere retroactivamente el supervisor nuevo.

## Carpetas cortas y progreso

El batch reserva `output/batch_1/` o el siguiente correlativo. Dentro usa identificadores por entrada e intento, como `g1_1`, conservando el nombre completo en metadatos. El modo simple, fuera del batch, sanea y trunca el identificador base a 12 caracteres antes de agregar el correlativo.

Durante el procesamiento se imprime `Analizando gen 4/20: NombreGen`. Los detalles y advertencias individuales se redirigen a logs, evitando veinte reportes en consola. La salida final muestra cantidades, ranking compacto y carpeta de resultados, con la aclaración heurística.

## Mapa de archivos

Ejemplo estructural de un batch; las cantidades, referencias y números de intento varían:

```text
output/batch_1/
├── resultados_genes.csv
├── ranking_genes.csv
├── resumen_batch.json
├── g1_1.log
├── inputs/g1_1.fasta
└── output/g1_1/
    ├── candidate.json
    ├── designer_handoff/
    ├── inputs/
    ├── pipeline_result.json
    ├── reference_comparison.json
    ├── all_references_windows.csv
    └── references/01_GITV.../
        ├── candidate.json
        ├── designer_handoff/run.json
        ├── inputs/designer_original.json
        ├── inputs/designer_original.json.bowtie.tmp
        ├── inputs/green_mask.json
        ├── offtarget_runs/
        ├── offtarget_selected_window/
        ├── offtarget_w2/             # Solo si se prueba otra candidata
        ├── integrated/
        ├── window_selection/
        ├── performance.json
        └── pipeline_result.json
```

| Archivo o carpeta | Contenido |
|---|---|
| `candidate.json` | Secuencia normalizada, nombre original y metadatos |
| `designer_handoff/run.json` | Parámetros, runtime y referencias a entradas del motor |
| `inputs/designer_original.json` y evidencia Bowtie | Registros si-Fi y alineamientos conservados |
| `inputs/green_mask.json` | Posiciones verdes, contexto de main target y advertencias |
| `integrated/report.json` | Candidato, métricas originales, evaluaciones, procedencia y advertencias |
| `integrated/records.csv`, `positions.csv`, `evaluations.csv`, `offtarget_matches.csv` | Tablas de registros, observaciones y alineamientos |
| `window_selection/all_windows.csv` | Todas las ventanas y métricas escalares, no solo la ganadora |
| `window_selection/selection_report.json` | Método, estado, procedencia e intentos del filtro |
| `window_selection/selected_window.json` | Solo existe si se exportó una ganadora |
| Figuras `.svg`/`.png` | Visualizaciones complementarias; no sustituyen la evidencia tabular |
| `offtarget_runs/` | Corridas de humano, maíz y abeja sobre el gen completo |
| `offtarget_selected_window/`, `offtarget_wN/` | Evidencia por intento de candidata |
| `reference_comparison.json`, `all_references_windows.csv` | Consolidación multirreferencia del gen |
| `performance.json`, `headless_performance.json` cuando existen | Tiempos por etapas y estado de caché; no una validación científica |

Por compatibilidad, `offtarget_selected_window/` corresponde a la **primera candidata probada**, aunque sea rechazada. La carpeta de la aprobada se identifica en `offtarget_selection.attempts`; no se debe deducir la ganadora del nombre de esa carpeta.

`pipeline_result.json` tiene distinto contenido según su nivel: en la raíz multirreferencia es una consolidación; dentro de cada referencia resume su pipeline. Revise campos de esquema y contexto.

## Tablas generales

`resultados_genes.csv` conserva una fila por entrada, en el orden original. Incluye estado y error, referencia destacada, coordenadas, longitud, clase, puntaje original, penalización, ajustado, calidad de evidencia, GC, cantidades de siRNAs, porcentajes de eficientes/diversidad y estados off-target.

`ranking_genes.csv` ordena primero los resultados evaluables por puntaje ajustado descendente; en empate conserva el orden de entrada. Los no evaluables quedan al final como `no_evaluable_result`, sin rango ni puntaje inventado. Un **cero calculado** sigue siendo evaluable; un nulo no se reemplaza por cero.

La selección para ranking requiere aprobación del filtro de la ventana. Puede haber un resultado de gen `partial` con una ventana aprobada entre las referencias completadas. El batch registra ese gen como `failed` si el workflow no completó sin errores, pero puede conservar su resultado evaluable en el ranking; el detalle `partial` permanece en el JSON individual y en `resumen_batch.json`. No confundir `complete` del proceso con eficacia biológica, ni asumir que una ejecución completada siempre tenga ganadora.

Las clases corta/convencional, suficiencia de evidencia, estados por especie y liberación de secuencia permanecen en columnas separadas. El ranking general es una comparación heurística entre contextos; no es una predicción validada de eficacia ni prueba que los genes sean objetivos equivalentes.

## Errores y reanudación

```powershell
python -m rnai_designer batch --resume output/batch_1
```

La reanudación utiliza las entradas normalizadas guardadas en `resumen_batch.json`, sin requerir el FASTA original. Opcionalmente:

```powershell
python -m rnai_designer batch genes.fasta --resume output/batch_1
```

En ese caso se verifica que el hash del FASTA coincida. No cambiar el archivo para quitar/agregar genes a un batch existente: iniciar otro batch con el archivo modificado.

Los genes completados con reporte íntegro no se repiten. Si el reporte final se escribió antes de una interrupción del registro batch, puede recuperarse. Los fallidos/interrumpidos se ejecutan en una carpeta nueva, conservando intentos anteriores; no se garantiza reanudación a mitad de un GITV. Los inválidos permanecen registrados y requieren un nuevo batch para corregir la entrada.

Existe bloqueo del sistema para impedir dos procesos sobre el mismo batch. Los CSV y el resumen se actualizan atómicamente después de cada gen; no existe una transacción única sobre todos los archivos. Código de salida del batch: **0** sin errores, **2** con errores y **130** al capturar Ctrl+C.

Las políticas y el código pueden haber cambiado desde un reporte antiguo. La reanudación no recalcula automáticamente cada gen marcado como completo ni reescribe históricamente todas las decisiones anteriores; documente con qué revisión se generó cada intento.
