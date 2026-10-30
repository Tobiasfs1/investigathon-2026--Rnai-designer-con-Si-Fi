# Interpretación y criterios científicos

[Índice](README.md) · [Arquitectura](arquitectura.md) · [Reportes](batch-y-reportes.md)

## Referencias GITV y colores

Un identificador como `GITV01000968_1` es el nombre de una secuencia representada en el índice Bowtie consultado. El programa no demuestra por ese nombre que se trate de un gen único, una isoforma concreta o una referencia taxonómicamente correcta. Un gen de entrada puede coincidir con numerosos GITV, por homología, redundancia u otras características del índice; la interpretación requiere su FASTA y anotación de origen.

Para cada main target, si-Fi distingue sus hits de los hits a otras referencias. Las **zonas verdes** se reconstruyen con la máscara de `general_helpers.get_target_data` y la lógica original de dibujo. Las zonas **rojas/off-target** corresponden al contexto de alineamientos a referencias no elegidas dentro de esa base primaria. No son equivalentes al filtro independiente contra humano, maíz y abeja.

La reconstrucción es sensible al orden de las filas: si-Fi acumula posiciones y guarda snapshots de main target restando los off-target ya vistos. No es correcto sustituirla por una simple resta global de conjuntos ni reorganizar los hits. También puede acortar o suprimir una región cuando una posición tiene hits tanto main-target como off-target.

**Verde no significa eficacia experimental ni seguridad general.** Las etiquetas dependen de los datos presentes y de la referencia elegida.

## Coordenadas, observaciones y ventanas

- `start_1based` y `end_1based_inclusive` expresan coordenadas inclusivas sobre la secuencia query del candidato.
- No son coordenadas genómicas ni posiciones del transcrito GITV. Bowtie conserva por separado offsets de referencia de base 0.
- Cada análisis GITV conserva su propio contexto y máscara; no mezclar sus anotaciones por compartir coordenadas numéricas.
- Un siRNA participa en una ventana cuando su intervalo completo está contenido en ella; no se interpola evidencia ausente.
- Las observaciones se agrupan por `(sirna_name, sirna_position, sirna_sequence)`. Varios hits del mismo siRNA no deben contarse como mediciones independientes de eficiencia.
- Si una métrica discrepa entre hits, se marca conflicto y se excluye ese valor de los componentes evaluables; los faltantes no se convierten en cero.

Configuración predeterminada: todas las longitudes enteras **180–250 nt**, con paso de **1 nt**, contenidas por completo en cada región verde. Una región verde menor que la longitud mínima configurada se conserva como una ventana exacta, en vez de alargarla artificialmente.

| Clase | Significado |
|---|---|
| `conventional_180_250nt` | Ventana dentro de 180–250 nt |
| `short_green_region` | Región verde completa menor que el mínimo configurado |
| `configured_range` | Ventana del intervalo configurado que no está en 180–250 nt |

Con valores predeterminados, una región verde de 27 nt pertenece a `short_green_region`. Cambiar el mínimo puede cambiar la clasificación; la penalización conserva su referencia fija de 180 nt.

## Puntaje original: 75 % eficiencia y 25 % GC

El código de [selección](../rnai-designer/rnai_designer/automatic_window_selection.py) define:

```text
E = siRNAs con is_efficient=True / siRNAs con clasificación evaluable
G = siRNAs con GC entre 30% y 52% / siRNAs con GC válido
D = secuencias siRNA distintas / observaciones siRNA agrupadas

efficacy_priority_score = 100 × (0.75 × E + 0.25 × G)
```

**E aporta 75 % y G aporta 25 %** del puntaje original. **D tiene peso 0 %**: conserva su cálculo, exportación y papel en los desempates. Se mantienen las comprobaciones existentes de disponibilidad de métricas; no se imputan valores faltantes. Los denominadores pueden diferir: clasificaciones faltantes no entran en el denominador de `E`.

El GC mostrado para la **ventana** es un dato descriptivo distinto de `G`, que cuenta **siRNAs** dentro del intervalo GC objetivo. En el cálculo de GC de siRNAs se consideran bases canónicas A/C/G/T; si no existen, el valor no es evaluable. La entrada acepta IUPAC, pero la presencia de bases ambiguas puede impedir acreditar el filtro off-target.

Los campos relevantes incluyen `is_efficient_true_count`, `is_efficient_evaluable_count`, `sirna_unique_count` y `sirna_unique_sequences_count`. La diversidad usa secuencias distintas sobre observaciones agrupadas, no sobre todas las filas crudas de hits.

## Penalización por longitud

Definición en [priority_scoring.py](../rnai-designer/rnai_designer/priority_scoring.py):

```python
length_penalty = weight * max(0, 1 - length_nt / 180)
adjusted_priority_score = max(0, efficacy_priority_score - length_penalty)
```

`weight` vale **30** por omisión y se configura en `workflow` mediante `--length-penalty-weight`. Debe ser finito y no negativo. La penalización resta puntos, no un porcentaje relativo del puntaje. Desde 180 nt vale cero. Si el puntaje original es nulo, el ajustado también permanece nulo.

Ejemplo aritmético: 27 nt con puntaje original 66,6667 recibe 25,5 puntos de penalización y queda en **41,1667**. Los nombres de campos y los componentes se conservan en los reportes. Los resultados anteriores y las auditorías históricas no se recalculan: conservan la fórmula utilizada en su ejecución.

## Calidad de evidencia y desempates

`quality_score` es una auditoría de documentación, no eficacia ni un cuarto componente del puntaje:

```text
100 × media(
  completitud de métricas,
  consistencia entre coordenadas y secuencia,
  cobertura de posiciones de inicio,
  integridad frente a duplicados exactos,
  trazabilidad
)
```

Los componentes se exportan por separado. La trazabilidad revisa disponibilidad de origen/versión de referencia, hash de la fuente, query ID y tamaño de siRNA registrado.

Se ordena por puntaje ajustado descendente. En empate se consideran, en este orden: fracción eficiente, fracción de GC objetivo, diversidad, calidad de evidencia, inicio menor, longitud menor y `region_id`. Para comparar ganadoras entre referencias se aplica además el identificador de referencia como desempate final. No se modifica el puntaje para representar esos desempates.

## Selección final y filtro off-target

La política `ranked_until_offtarget_pass_v1` recorre el ranking hasta encontrar una candidata que acredite el criterio. No basta con encabezar la prioridad. Un rechazo, error o evaluación indeterminada se registra y se continúa. Si ninguna pasa, `selected_window` queda nulo y no se publica secuencia ganadora.

Para cada ventana, el criterio requiere:

- Cero coincidencias exactas de **21 nt** en **Homo sapiens, Zea mays y Apis mellifera**.
- Tres evaluaciones con evidencia suficiente y resultados verificables.
- Bowtie exhaustivo, cero mismatches, sin gaps, todas las alineaciones y ambas hebras; parámetros y finalización documentados.
- Correspondencia de query/hashes, procesamiento de todas sus ventanas de 21 nt y archivos de resultados consistentes.

La ventana menor de 21 nt no puede acreditar este filtro. El filtro del **gen completo** puede detectar coincidencias fuera de una ventana que sí cumpla el criterio; deben revisarse los alcances por separado.

No se resta puntaje por fallar: los valores permanecen disponibles y cambia la elegibilidad de la ganadora. Los reportes registran `selected_priority_rank`, intentos y motivos.

### Estados que no deben confundirse

| Campo/estado | Interpretación |
|---|---|
| Criterio `criterion_met` | Cumple el criterio computacional documentado |
| Criterio `exact_match_detected` | Hay al menos una coincidencia exacta |
| Criterio `indeterminate` | La evidencia no alcanza para acreditar el criterio |
| `not_evaluated` | No se realizó o no consta la evaluación |
| `error` | Fallo de ejecución o procesamiento |
| Estado de corrida `imported_unvalidated` | Etiqueta de procedencia/importación; revisar el criterio independiente |
| Filtro `passes: true` | Autorización final de la ventana, junto con `status: criterion_met` |

El campo de corrida `imported_unvalidated` puede coexistir con un criterio posteriormente acreditado: no debe interpretarse aisladamente. El código de evaluación distingue procedencia de la importación, evidencia de la corrida y aprobación del filtro.

## Límites de interpretación

Este criterio solo considera los índices y parámetros registrados. No descarta off-targets por mecanismos distintos, coincidencias parciales, efectos dependientes de dosis, otros organismos o secuencias ausentes del índice. Las etiquetas de eficiencia de si-Fi y el promedio heurístico no demuestran silenciamiento in vivo. Regiones cortas y ventanas convencionales siguen siendo clases distintas aunque se comparen en el mismo ranking.
