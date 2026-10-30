# Rendimiento, almacenamiento y problemas conocidos

[Índice](README.md) · [Arquitectura](arquitectura.md) · [Batch y reportes](batch-y-reportes.md)

## Reutilización implementada

| Optimización | Alcance | Qué se conserva separado |
|---|---|---|
| Sesión temporal si-Fi por gen | Enumeración y referencias comparten alineamientos Bowtie; RNAplfold y registros base se calculan cuando hacen falta y se reutilizan | Etiquetas `is_off_target`, máscara original, selección y archivos de cada GITV |
| Memoización de eficiencia | Argumentos idénticos dentro de una configuración fija del pipeline | Cada fila original y sus métricas |
| Caché off-target de sesión | Corrida del gen completo retenida durante los GITV; hasta cuatro corridas adicionales de ventanas | Evidencia copiada y estado de cada evaluación |
| Caché opcional `workflow --cache-dir` | Artefactos headless identificados por contenido, parámetros y dependencias | Clave incluye main target; no confundir con la sesión compartida |
| Índices y tablas compactas | Agregaciones de ventanas y columnas escalares | Fórmulas y filas exportadas |
| Ranking diferido | Prueba primero la mejor ventana conocida; ordena las restantes solo si falla | Orden y criterios de desempate |
| Reutilización de parsing | El JSON ya validado y la preparación de siRNAs no se cargan/calculan dos veces en esas etapas | Validación de contenido y procedencia |
| Exportación incremental | Consolidación multirreferencia mediante archivo temporal, evitando acumular todas las ventanas de todos los GITV en RAM | CSV completo |
| Referencias sin verde | Evita índices y cálculos exclusivos de ventanas inexistentes | Reportes, advertencias y evaluación del gen completo |

La sesión compartida es privada, secuencial y temporal. Se comprueban secuencia, ID, parámetros, rutas y firmas de dependencias; cambios durante la ejecución invalidan la reutilización. No se comparte entre genes ni constituye una caché persistente recuperable después de terminar el proceso. Una finalización abrupta puede dejar temporales del sistema.

La máscara continúa usando el helper original y el orden de sus filas. No se sustituyen los cálculos científicos ni se eliminan referencias por número de hits. El código original de terceros se mantiene separado del adaptador.

## Costos que siguen existiendo

- Arranque de Conda/Python legado por referencia, aunque se compartan cálculos.
- Validación, etiquetado y exportación de datos por cada GITV.
- JSON, CSV, gráficos y copias de evidencias individuales, además de consolidaciones.
- Búsquedas adicionales para ventanas que todavía no tengan un resultado reutilizable y verificable.
- Serialización del estado batch y de los intentos; muchos descartes aumentan el detalle almacenado.

La caché reduce ejecuciones de herramientas, **no elimina las copias de evidencia en disco**. No hay porcentajes de aceleración o tamaños máximos medidos que puedan generalizarse a cualquier gen.

## Crecimiento con muchas referencias GITV

Para una región verde de longitud `L`, el número de ventanas convencionales con los parámetros predeterminados es:

```text
sum(max(0, L - longitud + 1) para longitud = 180..250)
```

Se agrega una ventana exacta por región más corta que el mínimo. La suma se repite por región y por referencia. Las ventanas se solapan y no son observaciones experimentales independientes.

Con cientos de GITV, cada referencia puede conservar el conjunto grande de registros primarios, evidencias, tablas y figuras. Los consolidados vuelven a representar información a otros niveles; `records.csv` incluye además serialización de registros. Esto puede multiplicar el almacenamiento aunque Bowtie y RNAplfold ya no se recalculen por referencia.

El recorrido off-target **no tiene un límite artificial de intentos**. Si muchas candidatas fallan o no pueden verificarse, se pueden crear numerosas carpetas `offtarget_wN/`, alineamientos, manifiestos e intentos JSON. Las reejecuciones fallidas del batch conservan sus carpetas antiguas. La nueva no reemplaza a la anterior.

Estos son problemas actuales, no funciones resueltas por la caché. Una futura deduplicación física de artefactos o un formato de almacenamiento diferente requeriría diseño, compatibilidad y validación específica; no está garantizada ni activada por esta documentación.

## Restricciones de Windows

El comando simple usa un identificador base saneado de hasta 12 caracteres más correlativo; batch utiliza identificadores cortos por entrada e intento. Se comprueba la longitud de ruta con margen para archivos internos. Esto reduce riesgos, pero no garantiza compatibilidad universal: la profundidad del directorio del proyecto, nombres GITV y rutas temporales también influyen. El código legado y las DLL pueden tener restricciones adicionales.

Utilice una ubicación corta al preparar una nueva ejecución y compruebe espacio libre suficiente. No mover, renombrar, comprimir ni borrar carpetas que estén siendo escritas por un batch activo. La limpieza de resultados requiere una decisión explícita sobre qué evidencia conservar y verificar primero el estado del proceso.

## Otros problemas y límites conocidos

- Python 2.7 y PyQt4 están obsoletos y no tienen una receta reproducible completa en esta copia.
- La copia original de `data_to_json` exporta los hits de hebra positiva; no basta para concluir ausencia de off-target. El filtro independiente consulta ambas hebras.
- Las ramas originales con ausencia de hits requieren cautela; vacío o incompleto no significa resultado negativo.
- Versiones y origen de las bases y RNAplfold no están totalmente documentados; no inferirlos de fechas o nombres.
- Una carpeta con archivos o un estado de proceso completado no implica ganadora aprobada.
- Los reportes y documentación de revisiones previas pueden describir selección de una sola candidata; no asumir que fueron migrados.
- No se ha demostrado eficacia biológica ni seguridad general con tests unitarios o una corrida sintética.

## Pruebas y medición

Las comprobaciones existentes incluyen datos pequeños, comparación de registros/máscaras con métodos originales, invariancia de métricas y pruebas de número de llamadas con motores simulados. Eso aporta evidencia de equivalencia local, no una medición de aceleración sobre 224 GITV.

Durante esta revisión documental no se ejecutaron ni tests ni benchmarks. No se lanzó el batch real. Para evaluaciones posteriores, separar tiempo de herramientas, parsing, exportación y tamaño en disco, registrar hardware/entorno y evitar competir con otra ejecución activa.
