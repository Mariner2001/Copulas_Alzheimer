# Resultados agregados

Salidas del pipeline que **sí** pueden publicarse: son estadísticos agregados, no datos
individuales. Se incluyen para que quien lea el artículo pueda verificar las cifras sin
solicitar acceso a ADNI.

| Fichero | Contenido |
|---|---|
| `bootstrap_B1000.csv` | Frecuencias de inclusión, errores de Monte Carlo, intervalos de τ y κ, familias, por par |
| `tail_evidence.csv` | ΔBIC, Vuong, cuotas de esquina, referencia gaussiana, cocientes |
| `strat_evidence.csv` | Los mismos contrastes dentro de cada estratificación |
| `stratified.csv` | Familia, τ y λ_L por análisis y par |
| `jitter_agg.csv` | Estabilidad solo-jitter frente a bootstrap |
| `gcgm_results.csv` | GCGM: EIP y correlación parcial de los 91 pares |
| `gcgm_summary.txt` | Resumen del GCGM con las cifras citadas en el artículo |
| `flow.csv` | Flujo de exclusión de la muestra (5146 → 630) |
