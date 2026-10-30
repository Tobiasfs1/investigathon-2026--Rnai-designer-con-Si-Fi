# Auditoría técnica del batch 1

[Guía del proyecto](../README.md) · [Detalle JSON](batch_1_auditoria.json) · [Tabla de 49 referencias](batch_1_referencias.csv)

## Alcance y conclusión

Auditoría de lectura sobre los ocho genes existentes en `rnai-designer/output/batch_1`. No se ejecutaron Bowtie, RNAplfold, el workflow ni el batch. No se alteraron código, criterios, selección ni archivos de resultados.

**No se detectaron errores de coordenadas, truncamiento de regiones verdes exportadas, aplicación de penalización ni propagación individual→batch. Sí se confirmó una limitación de implementación del si-Fi original: su exportador descarta los alineamientos negativos. Por eso las ventanas cortas no prueban una carencia biológica de homólogos largos.**

Las máscaras coinciden exactamente con el helper original, incluido su comportamiento dependiente del orden. Esa compatibilidad no valida científicamente todos los supuestos del helper.

## Conteos por gen

Las regiones se cuentan por contexto de referencia, no como unión deduplicada entre GITV. El máximo es la longitud de un tramo continuo, no la distancia entre el primer y último nucleótido verde.

| Gen | GITV completos | GITV con verde | Regiones continuas | Máximo verde nt | Ventana seleccionada nt | Ventanas convencionales evaluadas |
|---|---:|---:|---:|---:|---:|---:|
| Rpt5 | 2 | 2 | 2 | 34 | 34 | 0 |
| EloC | 15 | 5 | 9 | 409 | 194 | 26309 |
| mts | 3 | 1 | 1 | 21 | 21 | 0 |
| Hsp83 | 10 | 1 | 1 | 27 | 27 | 0 |
| snRNP-U1-70K | 5 | 2 | 2 | 22 | 22 | 0 |
| Rpt4/Rpt4R | 3 | 3 | 11 | 68 | 22 | 0 |
| Rpt3/Rpt3R | 6 | 1 | 2 | 1107 | 183 | 72562 |
| VhaSFD | 5 | 2 | 7 | 176 | 38 | 0 |

Total: **49 GITV, 35 regiones verdes y 98.901 ventanas revisadas**, de las cuales **98.871 son convencionales**. Las 49 referencias figuran completas. Solo EloC y Rpt3 tienen verde continuo >=180 nt.

## Verificaciones realizadas

- Ejecución aislada del cuerpo original de `general_helpers.get_target_data` con los JSON guardados: se adaptó únicamente `range` a la semántica de lista de Python 2. No se importó la GUI ni el pipeline científico.
- Comparación posición por posición entre verde original, sidecar, reconstrucción optimizada, regiones de selección y resumen por referencia. El rojo se contrastó con la unión original de intervalos off-target.
- Comparación ordenada de cada registro exportado con los alineamientos Bowtie positivos guardados. Los negativos permanecen en `.bowtie.tmp` y no en el JSON.
- Comprobación de coordenadas y secuencias de todos los registros y ganadoras; pertenencia de ganadoras a verde; recálculo de componentes y conteos de las ventanas elegidas desde registros originales.
- Verificación de la enumeración completa de cada CSV según región, rango y paso: 98.901 filas, sin faltantes, duplicados estructurales ni regiones omitidas. Comparación celda a celda de los CSV consolidados con sus fuentes individuales.
- Comprobación aritmética de los puntajes y la penalización en todas las filas, más comparación de ganadoras individuales, consolidados, snapshot batch, resultados y ranking.
- Verificación del hash registrado del resultado por gen y de identidad de alineamientos y métricas base entre referencias; solo cambian las etiquetas dependientes del target.

No se hallaron inconsistencias de estos tipos. Los resúmenes añaden deliberadamente `score_name` y `score_interpretation` frente a `selected_window.json`; los campos compartidos son idénticos.

## Diagnóstico por gen

### Rpt5

Entrada `cand0002_Rpt5 tcas=656024 StringTie STRG.22351.1 (coincide con Helixer) intrones=7/7 cds_cubierto=1.0`. Ganadora `GITV01002308_1`, 1417–1450, 34 nt.

**Confirmado:** Las dos máscaras originales suman dos tramos (4 y 34 nt); no existe región verde de 180 nt en estos resultados. La ganadora de 34 nt coincide con su secuencia y puntaje.

**Requiere revisión:** Hay cobertura main positiva continua de 705 nt en GITV01000966_1, reducida por hits a otras referencias bajo la política de un main target por vez. Revisar significado biológico de esas referencias y alcance de exclusividad, no asumir ausencia de homología larga.

### EloC

Entrada `cand0005_EloC tcas=664284 StringTie STRG.15183.1 (coincide con Helixer) intrones=2/2 cds_cubierto=1.0`. Ganadora `GITV01017907_1`, 1608–1801, 194 nt.

**Confirmado:** Regiones verdes máximas de 409, 378 y 219 nt; las 26.309 ventanas convencionales están presentes. La ganadora de 194 nt es consistente.

**Requiere revisión:** Un artefacto original deja la posición 868 simultáneamente verde y roja en GITV01000823_1. No afecta a la ganadora; revisar dependencia del orden de filas si se cambia la implementación.

### mts

Entrada `cand0007_mts tcas=662352 StringTie STRG.7615.1 (coincide con Helixer) intrones=5/5 cds_cubierto=1.0`. Ganadora `GITV01004263_1`, 428–448, 21 nt.

**Confirmado:** Solo un tramo verde exportado de 21 nt y un siRNA evaluable, marcado no eficiente. No hay ventanas convencionales omitidas de esa máscara.

**Requiere revisión:** GITV01003170_1 tiene 1.135 alineamientos negativos y cobertura continua de 315 nt (974–1288), no exportada por data_to_json. Revisar limitación de hebra y orientación de las referencias/candidato antes de concluir que no hay regiones largas.

### Hsp83

Entrada `cand0010_Hsp83 tcas=656270 StringTie STRG.14232.1 (coincide con Helixer) intrones=2/2 cds_cubierto=1.0`. Ganadora `GITV01022783_1`, 410–436, 27 nt.

**Confirmado:** Solo un tramo verde exportado de 27 nt; los siete siRNAs evaluables tienen is_efficient=False. Penalización 25,5 y ajustado 41,1667 correctos.

**Requiere revisión:** GITV01000126_1 contiene 2.404 alineamientos negativos y cobertura continua de 1.107 nt (1398–2504), ausente del JSON original por el filtro de hebra positiva. Revisión de orientación/exportación prioritaria.

### snRNP-U1-70K

Entrada `cand0011_snRNP-U1-70K tcas=662320 StringTie STRG.15363.1 (coincide con Helixer) intrones=8/8 cds_cubierto=1.0`. Ganadora `GITV01037515_1`, 575–596, 22 nt.

**Confirmado:** Dos tramos verdes de 21 y 22 nt; la ganadora de 22 nt tiene dos siRNAs evaluables, ninguno eficiente. Coordenadas, secuencia y puntaje consistentes.

**Requiere revisión:** GITV01001014_1 contiene 1.258 alineamientos negativos; el tramo continuo 1–390 no llega al JSON original. Revisar la política de hebras. No hay evidencia de que solo existan coincidencias cortas.

### Rpt4/Rpt4R

Entrada `cand0015_Rpt4,Rpt4R tcas=656146 StringTie STRG.8886.1 (coincide con Helixer) intrones=8/8 cds_cubierto=1.0`. Ganadora `GITV01001584_1`, 855–876, 22 nt.

**Confirmado:** Los 11 tramos verdes están presentes, incluido 1293–1360 (68 nt). La ventana 855–876 (22 nt) obtiene 23,6667 ajustados frente a 20,2222 de la de 68 nt: no fue omitida. La ganadora tiene 0/2 siRNAs eficientes.

**Requiere revisión:** La cobertura main positiva llega a 1.360 nt, pero otras referencias restan casi toda ella con la lógica original. Además hay dos posiciones verdes/rojas aisladas en otros GITV; ninguna cae en la ganadora. Revisar exclusividad entre referencias y adecuación biológica del puntaje.

### Rpt3/Rpt3R

Entrada `cand0021_Rpt3,Rpt3R tcas=656344 StringTie STRG.19698.1 (coincide con Helixer) intrones=8/8 cds_cubierto=1.0`. Ganadora `GITV01000607_1`, 321–503, 183 nt.

**Confirmado:** Las dos regiones verdes, de 1.107 y 343 nt, generan 72.562 ventanas convencionales y todas están exportadas. La ganadora de 183 nt es consistente y fue la primera candidata aprobada.

**Requiere revisión:** Hay también alineamientos negativos omitidos, aunque la hebra positiva ya permite regiones largas. Validar experimentalmente el resultado; el puntaje no es eficacia medida.

### VhaSFD

Entrada `cand0024_VhaSFD tcas=657235 StringTie STRG.6645.2 (coincide con Helixer) intrones=9/9 cds_cubierto=1.0`. Ganadora `GITV01006946_1`, 334–371, 38 nt.

**Confirmado:** La región de 176 nt está presente y también fue aprobada como ganadora de otra referencia. Obtiene 47,8376 ajustados, por debajo de los 48,5556 de la ventana de 38 nt. No hay regiones verdes >=180 nt exportadas.

**Requiere revisión:** Existen tramos negativos de 221 nt en GITV01003434_1 que no se exportan. Revisar hebras/orientación y el equilibrio del puntaje: la penalización se aplica pero no impone una longitud mínima ni favorece siempre a la más larga.

## ¿Hay regiones largas que faltan?

Sí existen coberturas de alineamientos main más largas que las máscaras, pero hay dos causas distintas:

1. **Exclusividad entre GITV:** en Rpt5 hay cobertura positiva continua de 705 nt y en Rpt4 de 1.360 nt. Otros hits primarios hacen que el helper reste casi toda esa cobertura al analizar cada GITV como main target separado. No son regiones verdes perdidas por el reportador; requieren interpretar biológicamente qué representan las referencias que se restan.
2. **Exportación unilateral heredada:** mts, Hsp83 y snRNP-U1-70K tienen referencias con cobertura larga exclusivamente negativa, omitida por `data_to_json`. VhaSFD también tiene tramos negativos de 221 nt. Es una limitación real del software original, no evidencia de que los datos carezcan de tramos largos.

Estos tramos de alineamientos no pueden declararse automáticamente específicos, eficientes ni aprobados para las tres especies. Revisar su orientación y la política de hebras exigiría una decisión científica y validación posterior; no se modificó nada aquí.

EloC conserva sus tramos de 409/378/219 nt y Rpt3 los de 1.107/343 nt. Los 68 nt de Rpt4 y los 176 nt de VhaSFD también están en las tablas. Ninguno se perdió por fragmentación accidental en la exportación.

## Verde/rojo: artefactos originales

Se detectaron tres posiciones verdes que también pertenecen al rojo final acumulado: EloC/GITV01000823_1: **868**; Rpt4/GITV01003734_1: **1292**; Rpt4/GITV01016098_1: **1037**. Es la consecuencia exacta de guardar snapshots antes de conocer todos los off-target posteriores. La reconstrucción de rnai-designer reproduce esas posiciones, no las inventa. Son tramos aislados, no regiones largas ocultas, y ninguna ventana ganadora los contiene.

El script auxiliar `siFi21--1.2.3/region_verde.py`, que no es utilizado por este workflow, además usa min/max para unir potenciales huecos y aplica un desplazamiento de base de coordenadas. No utilizar una salida de ese script como prueba de una región verde continua sin revisar estos puntos.

## Penalización y elección de longitud

`length_penalty = 30 * max(0, 1 - length_nt/180)` y `adjusted = max(0, original - penalty)` se cumplen en las 98.901 filas y en los dos niveles individual/batch. La comparación usa tolerancia numérica absoluta 1e-9 y relativa 1e-11 para las representaciones de punto flotante.

| Gen | Longitud | Original | Penalización | Ajustado | Eficientes/evaluables |
|---|---:|---:|---:|---:|---:|
| Rpt5 | 34 | 71.4286 | 24.3333 | 47.0952 | 3/14 |
| EloC | 194 | 56.8966 | 0.0000 | 56.8966 | 27/174 |
| mts | 21 | 33.3333 | 26.5000 | 6.8333 | 0/1 |
| Hsp83 | 27 | 66.6667 | 25.5000 | 41.1667 | 0/7 |
| snRNP-U1-70K | 22 | 66.6667 | 26.3333 | 40.3333 | 0/2 |
| Rpt4/Rpt4R | 22 | 50.0000 | 26.3333 | 23.6667 | 0/2 |
| Rpt3/Rpt3R | 183 | 64.0082 | 0.0000 | 64.0082 | 16/163 |
| VhaSFD | 38 | 72.2222 | 23.6667 | 48.5556 | 3/18 |

La penalización es correcta, pero no impone longitud mínima. En Rpt4, 22 nt obtiene **23,6667** frente a **20,2222** de 68 nt. En VhaSFD, 38 nt obtiene **48,5556** frente a **47,8376** de 176 nt. Es el resultado de la fórmula existente, no un fallo de aplicación. En cuatro genes la ganadora contiene cero siRNAs eficientes; GC y diversidad mantienen un puntaje positivo porque no existe un umbral mínimo de eficiencia.

Los **11 intentos off-target por referencia con candidata** aprobaron su primera ventana. Ninguna selección corta de este batch fue causada por descartar primero una ventana larga en el recorrido off-target. Esto confirma el comportamiento registrado; no constituye una reejecución ni una nueva certificación biológica de las bases/herramientas.

## Software frente a limitaciones bioinformáticas

| Categoría | Resultado de auditoría |
|---|---|
| Implementación actual de máscaras, coordenadas y enumeración | Compatible con los artefactos y el helper original; sin discrepancias detectadas |
| Propagación de puntajes individual/batch | Correcta, sin omisión de penalización |
| Exportador original de si-Fi | Omite la hebra negativa; revisión prioritaria para mts, Hsp83 y snRNP-U1-70K, también relevante para VhaSFD |
| Helper original de colores | Dependiente del orden; tres bases con verde/rojo superpuesto, fuera de ganadoras |
| Especificidad respecto de otros GITV | Reduce tramos positivos largos en Rpt5 y Rpt4; revisar anotación y objetivo de exclusividad |
| Preferencia de la heurística | Puede favorecer regiones cortas y candidatos con cero eficientes; comportamiento especificado, adecuación biológica pendiente |
| Script auxiliar region_verde.py | Riesgo de unir discontinuidades y desplazar coordenadas; no intervino en el batch auditado |

## Pruebas pequeñas y límites

**8 pruebas unitarias aprobadas en 1,50 s**: reconstrucción ordenada y casos de longitud 27/100/179/180/200/250 nt, además de ausencia de puntaje. No se ejecutaron pruebas largas ni el batch.

La auditoría verifica integridad interna y compatibilidad con la lógica original. No demuestra identidad taxonómica, anotación de los GITV, calidad del ensamblado, validez experimental del score, ni qué selección resultaría tras una futura modificación de hebras o criterios. Los reportes originales permanecieron intactos; estos archivos son documentación adicional de auditoría.
