# Hebras negativas en si-Fi21: origen del filtrado y efecto en la GUI

[Auditoría previa del batch 1](batch_1_auditoria.md) · [Verificación mecánica y datos](hebras_sifi21_verificacion.json)

## Dictamen

La exclusión **existe en el código oficial y en la GUI original**; no fue introducida por rnai-designer ni por su caché. Es un filtro explícito en la conversión de alineamientos a JSON. No ocurre en el helper que construye `main_target_dict`.

La evidencia es compatible con una restricción de orientación del software: el manual pide el mRNA objetivo como query de diseño y el programa advierte cuando query/base no están orientadas correctamente. Sin embargo, no se encontró una explicación de los autores que declare y justifique científicamente la omisión de todos los hits negativos. Por tanto, no se puede asegurar la intención completa ni descartar un defecto heredado en la implementación. “Filtrado explícito y preexistente” no equivale a “decisión científicamente validada”.

Esto precisa la auditoría anterior: llamar a la omisión una limitación heredada describe su efecto real sobre estos datos; no demuestra que sea un accidente del autor ni que agregar filas negativas directamente sea correcto.

## Comparación con la distribución oficial

Se consultó el repositorio `snowformatics/siFi21-`, no solo la copia local. La etiqueta **v1.2.3-0008** apunta a `65c197250f396783022f07907ca44cf5509f24e0`. Los cuerpos de `data_to_json` y `get_target_data` coinciden exactamente con los métodos locales, normalizando únicamente finales de línea al leerlos.

El historial accesible de `sifi_pipeline.py` contiene un único commit, `5d86eada37bbaf929fb5a658aef3c0cbc3597520`, descrito como Initial commit. No aporta un mensaje que explique por qué se adoptó ese filtro. El método también coincide con master consultado durante esta revisión.

Fuentes oficiales:

- [Pipeline de la etiqueta v1.2.3-0008](https://github.com/snowformatics/siFi21-/blob/65c197250f396783022f07907ca44cf5509f24e0/sifi_pipeline.py).
- [Helper original de regiones](https://github.com/snowformatics/siFi21-/blob/65c197250f396783022f07907ca44cf5509f24e0/general_helpers.py).
- [Código de la GUI/gráficos](https://github.com/snowformatics/siFi21-/blob/65c197250f396783022f07907ca44cf5509f24e0/show_plot.py).
- [Manual de los autores](https://github.com/snowformatics/siFi21-/blob/7ef3e82003883156ecfa2cc73c5d36e93e02f9c2/Quick_help.pdf), páginas 2, 6 y 7.
- [Artículo de los autores](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2019.01023/full), secciones Implementation y Off-Target Searching and Prediction.

## Dónde se pierden las filas

| Etapa | Tratamiento de los negativos |
|---|---|
| `run_bowtie` | Su comando no pide restringir a una sola hebra; la evidencia guardada contiene ambas |
| `bowtie_to_lst` | Convierte todas las líneas; no elimina filas por signo |
| `get_main_target` | Cuenta nombres sobre la lista Bowtie completa: un GITV con solo negativos puede aparecer en el selector |
| `data_to_json` | El cálculo y `json_lst.append` están dentro del bloque `if strand == '+'`; el negativo no se añade |
| `general_helpers.get_target_data` | No consulta `strand`; separa main/off por `is_off_target` sobre el JSON recibido |
| `DrawPlot.plot_design` | Obtiene los diccionarios del helper y dibuja sus posiciones; no vuelve a leer Bowtie ni recupera los negativos |

Así, una referencia puede aparecer en el diálogo de selección pero terminar sin regiones verdes. Cuando todos los hits exportables desaparecen, `run_pipeline` puede devolver el mensaje sobre orientación en lugar de mostrar un gráfico normal.

El filtro afecta también las filas que habrían sido off-target, no solo las del main target. Esta observación corresponde a la ruta original si-Fi. No debe trasladarse sin más al control independiente de rnai-designer para humano/maíz/abeja, que ejecuta sus propias búsquedas documentadas de ambas hebras.

## Qué sugiere intención y qué no demuestra

El manual distingue el uso de un trigger existente y el diseño a partir del mRNA objetivo; describe bases típicamente de cDNA/transcriptoma. El propio programa menciona orientación en su manejo de resultados vacíos. Es razonable inferir que las secuencias orientadas eran parte del uso previsto.

Pero el manual consultado no contiene una regla explícita de exclusión universal de Bowtie negativo. El código conserva una posible rama de tratamiento inverso/complementario comentada, sin ejecutarla. Eso muestra que el problema estaba presente en el desarrollo, pero no permite reconstruir la decisión del autor. El artículo distingue selección de la hebra guía por criterios termodinámicos y uso de bases transcriptómicas; no justifica por sí solo el filtro de orientación observado.

`strand` de Bowtie y la opción GUI “Strand selection” no son el mismo concepto. El primero describe la orientación del alineamiento. La opción interviene en las reglas de selección/eficiencia de la guía. Desactivarla no saca `json_lst.append` del bloque positivo.

## Datos reales conservados

Los siguientes son ejemplos de referencias individuales del batch existente; los conteos no son el total de hits del gen contra todas las referencias:

| Gen / referencia | Hits Bowtie + | Hits Bowtie − | Filas del main target en JSON | Verde de esa referencia |
|---|---:|---:|---:|---:|
| cand0007_mts / GITV01003170_1 | 0 | 1135 | 0 | 0 nt |
| cand0010_Hsp83 / GITV01000126_1 | 0 | 2404 | 0 | 0 nt |
| cand0011_snRNP-U1-70K / GITV01001014_1 | 0 | 1258 | 0 | 0 nt |

Los JSON de esas ejecuciones contienen respectivamente 1, 7 y 3 filas positivas de otras referencias. No están vacíos en conjunto, pero ninguno conserva una fila del target negativo de la tabla. Esto explica por qué sus regiones no aparecen como verdes al seleccionarlos. Los datos Bowtie originales permanecen disponibles.

## Verificación sin ejecutar el análisis

Se ejecutaron en memoria fragmentos puros del código oficial con entradas mínimas:

1. Dos filas Bowtie, una positiva y una negativa: `data_to_json` exporta solo la positiva, tanto con `strand_check=0` como con `strand_check=2`.
2. Una fila JSON negativa suministrada al helper: `get_target_data` construye sus 21 posiciones main. Demuestra que el helper no impone el filtro de hebra.
3. El selector de main targets recibe ambos nombres presentes en las filas Bowtie.
4. `plot_design` original, con llamadas gráficas y guardado interceptados, recibe el JSON filtrado y no dibuja verde para el target exclusivamente negativo.

Resultado: **5 comprobaciones mecánicas aprobadas**. La eficiencia se sustituyó por un retorno de prueba y la interfaz gráfica por objetos sin salida para aislar la selección de filas. No se calcularon ni validaron métricas científicas con esos sustitutos. No se abrió una ventana Qt ni se ejecutó el binario GUI instalado: la conclusión se refiere a la ruta comprobada del código fuente oficial y local.

El manual PDF se consultó en memoria; la documentación no se dedujo únicamente del mensaje de error. No se ejecutaron Bowtie, RNAplfold, el batch ni un workflow.

## Consecuencia para la interpretación

No atribuir automáticamente las ventanas cortas a falta de homología larga. Tampoco incorporar los negativos copiando las métricas positivas: la relación entre orientación, secuencia guía, posiciones y accesibilidad debe evaluarse antes. La revisión siguiente debe establecer cómo están orientados los candidatos y los transcritos GITV y si coinciden con el supuesto del motor.

Esta investigación no modifica código, clasificación de regiones, criterios ni resultados. Solo añade este informe y su evidencia JSON.
