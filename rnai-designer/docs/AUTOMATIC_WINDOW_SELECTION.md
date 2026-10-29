# Selector automático de ventanas

El comando `select-windows` recorre ventanas convencionales 180-250 nt con paso de 1 nt y agrega una ventana completa para cada run verde mas corto que 180 nt. Cada ventana permanece dentro de un unico run; no se extienden ni se unen regiones.

Las coordenadas exportadas son 1-based e inclusivas, igual que `sirna_position` y los runs reconstruidos. Internamente, al extraer la secuencia de Python, se convierten a slices 0-based/half-open. `workflow` conserva la secuencia seleccionada solo si supera el filtro descrito abajo. Ventanas repetidas se descartan por el par (inicio, fin); ventanas con longitudes distintas no son duplicadas.

## Regiones verdes cortas

Las regiones verdes menores de la longitud minima convencional se agregan como una sola ventana que cubre exactamente el run completo. Se etiquetan `short_green_region`, con `conventional_range_comparable=false`; las ventanas 180-250 nt se etiquetan `conventional_180_250nt`. No se extiende ni se une ningun run a traves de un hueco. Las columnas `sirna_record_count`, `data_sufficiency_status`, faltantes y clasificaciones indican si hay evidencia suficiente; una region corta sin siRNAs completos o sin `is_efficient` evaluable no recibe prioridad. Las coordenadas siguen siendo 1-based inclusivas.

## Métricas por ventana

`all_windows.csv` contiene una fila por ventana. Las columnas principales son:

| Métrica | Definición |
|---|---|
| `length_nt`, inicio y fin | Longitud exacta y coordenadas 1-based inclusivas |
| `gc_percent` | G+C sobre las bases canónicas A/C/G/T de la secuencia; bases ambiguas se excluyen del denominador; null si no hay bases canónicas |
| `sirna_record_count` | Filas originales de si-Fi cuyo intervalo completo cabe en la ventana, incluyendo filas repetidas/hits |
| `sirna_unique_count` | Grupos únicos por nombre, posición y secuencia, como en el explorador existente |
| `sirna_unique_sequences_count` | Secuencias distintas entre esos grupos; la misma secuencia en posiciones distintas se cuenta una vez aquí |
| `is_efficient_true_count` | SiRNAs etiquetados `true` por si-Fi |
| `is_efficient_evaluable_count` | Grupos con etiqueta original `true` o `false`; null y conflictos no son evaluables |
| `is_efficient_true_fraction_among_evaluable` | Conteo true dividido por evaluables; resumen de etiquetas importadas, no probabilidad |
| `sirna_complete_count` | SiRNAs con las cinco métricas originales disponibles y sin conflicto: `is_efficient`, `strand_selection`, `end_stability`, `target_site_accessibility`, `accessibility_value` |
| `missing_data_percent` | Celdas null o conflictivas / (grupos únicos × cinco métricas) × 100. Sin grupos, se informa 100% |
| Cobertura | Inicios únicos observados / inicios posibles de un siRNA completo; también se da la unión de nucleótidos cubierta por los intervalos y su fracción |
| `accessibility_value_*` | n, media, mediana, mínimo y máximo de los valores numéricos originales, excluyendo conflictos |
| Duplicados | Filas exactamente idénticas; filas adicionales debidas a múltiples hits; y siRNAs con conflictos de métricas |
| Off-target | Alineamientos cuyo intervalo de query se superpone con la ventana, con estado, trazabilidad y conteos separados por humano, maíz y abeja |

No se rellenan métricas ausentes ni se recalculan las etiquetas de si-Fi. El conteo de filas originales y el de siRNAs únicos están relacionados porque una secuencia puede aparecer en más de una fila por múltiples hits. El número de secuencias únicas no es el número de posiciones únicas. La completitud de métricas y el porcentaje faltante son complementarios; la fracción eficiente comparte el denominador de etiquetas evaluables. Cobertura de inicios y cobertura de nucleótidos también se relacionan, pero no son equivalentes: las ventanas siRNA pueden solaparse. Estas relaciones se mantienen visibles para que no se interpreten como evidencia independiente.

`accessibility_value` se copia sin transformar. El JSON integrado no conserva una unidad ni los detalles suficientes del parámetro de RNAplfold para asignarle aquí una interpretación física; por eso el reporte muestra el nombre original y estadísticas descriptivas, sin convertirlo a energía o probabilidad.

## Prioridad exploratoria y filtro de liberacion

La ventana se prioriza por un puntaje descriptivo separado del puntaje de calidad de evidencia. Para cada ventana se calculan tres fracciones en [0,1]:

```text
efficacy_priority_score = 100 / 3 * (
    si-Fi is_efficient == true / is_efficient evaluables
  + siRNAs con GC individual entre 30% y 52% / siRNAs con GC canonico disponible
  + secuencias siRNA distintas / posiciones siRNA unicas
)
```

Se requiere disponer de los tres componentes; no se imputan valores si-Fi faltantes. Los componentes quedan como columnas separadas en `all_windows.csv`. El GC por ventana completa tambien se conserva como estadistica descriptiva, pero el componente GC se calcula para cada siRNA de 21 nt contenido. El rango 30-52% corresponde a una regla clasica descrita por Reynolds et al. (Nature Biotechnology, 2004): https://doi.org/10.1038/nbt936. No esta validada aqui para Dalbulus maidis. La diversidad se define solo como proporcion de secuencias distintas entre siRNAs por posicion y no se interpreta como eficacia experimental.

`is_efficient` es la clasificacion original de si-Fi y ya combina las reglas que estaban activas para seleccion de hebra, estabilidad de extremos y accesibilidad. Por eso estas variables correlacionadas no vuelven a sumarse como predictores independientes. El promedio con pesos iguales es una regla transparente de priorizacion exploratoria, no un modelo biologico validado ni una prediccion de knockdown. El `quality_score` existente y sus componentes siguen disponibles para auditar documentacion, consistencia, cobertura y trazabilidad; no forman parte de la prioridad.

Empates: mayor fraccion `is_efficient`, mayor fraccion GC objetivo, mayor diversidad, mayor `quality_score`, inicio menor, longitud menor y `region_id`. Si las etiquetas `is_efficient` son faltantes/conflictivas en todas las candidatas, no se elige ni publica una ventana.

Antes de devolver la secuencia, `workflow` ejecuta Bowtie contra las tres bases configuradas usando los siRNAs de la ventana priorizada: longitud 21, cero mismatches, sin gaps, todas las alineaciones y ambas hebras. La secuencia se publica unicamente si cada ejecucion termina correctamente, acredita cobertura completa, el estado es `criterion_met` y hay cero hits exactos de 21 nt. Los reportes se guardan por especie en `offtarget_selected_window/`. Si aparece cualquier hit o una corrida es incompleta/indeterminada/error, se conservan coordenadas, metricas y auditoria, pero se omite la secuencia del resultado final. No se prueba automaticamente una ventana siguiente.

Este filtro significa ausencia computacional de coincidencias exactas de 21 nt en las secuencias indexadas; no es una certificacion general de seguridad biologica o ambiental. El orden puede ser distinto entre especies y la secuencia destacada requiere validacion experimental.
## Ejecutar y salidas

```powershell
python -m rnai_designer select-windows --report output/candidato/report.json --out output/candidato/selection
```

Se puede acotar el intervalo para revisión o pruebas con `--min-length` y `--max-length`; por omisión se usa 180–250 nt. El paso permanece fijado en 1 nt. Se acepta `--sirna-size` solo cuando el manifiesto no lo registra y no puede inferirse. La salida debe ser nueva o estar vacía.

- `selected_window.json`: coordenadas y auditoria; incluye la secuencia solo si pasa el filtro off-target estricto de las tres especies.
- `all_windows.csv`: tabla exhaustiva para comparación y auditoría.
- `selection_report.json`: rango, regiones y cantidad de ventanas, fórmula, desempate, procedencia, evaluación off-target por especie, advertencias y razón de selección.
- `selected_window.svg` y `.png`: gen completo, runs verdes, ventana priorizada en amarillo y tracks off-target por especie.

El selector no modifica archivos de si-Fi, bases de datos, la referencia ni las predicciones originales.
