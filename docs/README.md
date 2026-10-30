# Documentación técnica

[Inicio del proyecto](../README.md) · [CLI de rnai-designer](../rnai-designer/README.md)

## Índice

1. [Instalación y reproducibilidad](instalacion-y-reproducibilidad.md): versiones verificadas, entornos Python 3/Python 2, configuración y `doctor`.
2. [Arquitectura e integración](arquitectura.md): módulos, procesos y frontera con si-Fi original.
3. [Criterios científicos](criterios-cientificos.md): referencias GITV, colores, métricas, puntuación y off-target.
4. [Batch y reportes](batch-y-reportes.md): entradas, ranking, estados, archivos y reanudación.
5. [Rendimiento y limitaciones](rendimiento-y-limitaciones.md): optimizaciones actuales y crecimiento de almacenamiento.
6. [Créditos y licencias](creditos-y-licencias.md): componentes propios, terceros y permisos pendientes.

## Alcance de esta edición

Revisión documental del **29 de octubre de 2026**, basada en el código local de `rnai-designer` **0.3.0**, si-Fi **v1.2.3-0008**, dependencias, metadatos instalados y manifiestos de una ejecución histórica. No se ejecutaron el pipeline, `doctor`, herramientas científicas ni tests durante el batch activo.

La versión del paquete no identifica todas las revisiones locales: para reproducir una entrega, registre también el commit Git y los cambios sin commit. No se presupone que los reportes históricos hayan sido recalculados con la política actual de selección.

## Documentación especializada existente

Se conservan las ubicaciones originales:

- [Auditoría de si-Fi](../rnai-designer/docs/SIFI_ANALYSIS.md).
- [Integración](../rnai-designer/docs/SIFI_INTEGRATION.md) y [formatos detallados](../rnai-designer/docs/FORMATS.md).
- [Exploración de ventanas](../rnai-designer/docs/WINDOW_EXPLORATION.md) y [selección automática](../rnai-designer/docs/AUTOMATIC_WINDOW_SELECTION.md).
- [Automatización](../rnai-designer/docs/SIFI_AUTOMATION.md) y [portabilidad](../rnai-designer/docs/portabilidad.md).
- [Historial de optimización](../rnai-designer/docs/PERFORMANCE_OPTIMIZATION.md) y [validación](../rnai-designer/docs/VALIDATION.md).

Esos documentos incluyen decisiones y validaciones anteriores. Para selección sucesiva, caché por gen y batch actuales, utilice esta guía junto con el código. Las descripciones históricas no garantizan el comportamiento de revisiones posteriores.
