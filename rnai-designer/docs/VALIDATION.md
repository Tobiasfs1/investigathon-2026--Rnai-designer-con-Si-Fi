# Validación de la entrega

Fecha: 28/10/2026. Entorno: Windows, Python 3.14.10, matplotlib 3.11.2, pytest 9.2.0. Las dependencias completas resueltas para ese entorno están en `requirements-tested.txt`; los rangos soportados declarados están en pyproject.toml y requirements-dev.txt. No se ejecutó una matriz de versiones de Python.

Desde la raíz de rnai-designer:

```powershell
python -m pytest -q
python -m rnai_designer --candidate examples/candidate.json --sifi examples/sifi_result.json --run examples/run.json --offtargets examples/offtargets.json --out output/demo
```

Resultado actual: **100 pruebas aprobadas**, incluida una prueba de integracion headless con el runtime Conda Python 2 (`python -m pytest -q`), que cubren pruebas de reconstrucción compatible con la máscara original, sidecar explícito con varias regiones, workflow integrado, ventanas, datos faltantes/conflictivos, vínculo de intervalos off-target y exportaciones. También incluye un test que invoca Bowtie 0.11.3 contra índices sintéticos temporales; el test no modificó índices instalados. `rnai-designer doctor` detectó Bowtie 0.11.3, Bowtie-inspect y RNAplfold 2.1.5. `databases --inspect` verificó lectura y conteo de nombres en los cuatro índices configurados.

La búsqueda automática se probó con las bases sintéticas del test de integración, no con bases de organismos instaladas. No se ejecutó análisis del fixture artificial contra las bases biológicas. La CLI `analyze` debe usarse con secuencias reales del proyecto y una carpeta de salida nueva.

Las pruebas cubren conservación del candidato y prioridad, registros y extensiones originales, procedencia y parámetros, hits repetidos, métricas discordantes, null, eficiencia desactivada, coordenadas inválidas, alias de query, referencias discordantes, JSON inválido/no finito y claves repetidas, organismos duplicados, archivos ausentes y resultados vacíos. Las evaluaciones omitidas permanecen no evaluadas; los resultados vacíos nunca producen conclusión de seguridad.

Se compararon SHA-256 de los **75 archivos** de siFi21--1.2.3 antes y después del desarrollo inicial: sin modificaciones. El acceso de inspección de índices no escribe en ellos. No existe un repositorio Git en la raíz inspeccionada; por eso la comprobación de integridad se hizo mediante hashes.

Estas verificaciones demuestran funcionamiento del importador, el clasificador del criterio, el adaptador Bowtie y una ejecucion del motor original si-Fi en modo headless sobre un indice sintetico. **No son validacion biologica ni resultados de un gen candidato real.** Hay indices instalados para cuatro prefijos; faltan metadatos comprobados sobre release, accesiones de origen y exhaustividad taxonomica.

La nueva CLI `windows` se ejecutó además contra `output/demo/report.json`: reconstruyó una zona verde, generó tres ventanas de 21 nt y exportó JSON, CSV, SVG y PNG. La inspección visual confirmó que la ventana seleccionada queda marcada y que los ejes usan inicio 1-based. Son datos ficticios de ejemplo, no resultados biológicos.


La CLI select-windows recorre exhaustivamente las longitudes 180-250 nt en pasos de 1 nt, excluye ventanas que cruzan los limites de un run verde y exporta la tabla completa, la seleccion, el reporte metodologico y la visualizacion del gen. Los tests cubren ambos extremos, desplazamiento, gaps entre runs, ausencia de candidatas, empates, metricas faltantes/conflictivas, duplicados y separacion off-target por especie. La validacion automatizada usa fixtures sinteticos; no representa una seleccion biologica.

El comando `workflow` fue probado de extremo a extremo con el motor original y dos referencias Bowtie sinteticas: se generaron dos carpetas independientes, cada una con JSON original y mascara de get_target_data, y un reporte consolidado. Off-target se sustituyo por not_evaluated en esta prueba y no se consultaron indices de organismos.


La seleccion ahora usa una prioridad exploratoria con tres componentes auditables y pesos iguales: proporcion de etiquetas `is_efficient=True` originales entre las evaluables; proporcion de siRNAs individuales dentro del rango GC 30-52%; y proporcion de secuencias siRNA distintas entre posiciones. No se imputan etiquetas ausentes. La calidad de evidencia previa permanece independiente. Tests cubren calculo de componentes, faltantes y compuerta de liberacion.

`workflow` vuelve a ejecutar Bowtie para los siRNAs de la ventana priorizada, por especie, y solo conserva su secuencia en `selected_window.json` y `pipeline_result.json` si las tres evaluaciones completas registran cero hits exactos de 21 nt. Los tests de compuerta simulan estados aprobados/no evaluados; la integracion real del motor si-Fi conserva off-target en no evaluado y, correctamente, retiene la secuencia. No se ejecuto una secuencia candidata real contra los indices biologicos.


La ultima suite incluye enumeracion automatica de uno o dos main targets, preservacion separada de coordenadas, reporte consolidado y continuidad cuando falla una referencia; tambien cubre runs verdes menores de 180 nt, gaps y sidecars sin registros siRNA suficientes. La prueba de integracion invoco el motor si-Fi original con dos referencias sinteticas Bowtie desde un solo comando. En el fixture ambas referencias compartian las mismas guias: se verifico que la mascara del helper original puede reducirse/suprimirse por alineamientos off-target coincidentes y que el sidecar registra esta advertencia. Las metricas biologicas de ese fixture y sus filtros off-target no son resultados biologicos.
