# Resultados de E1 y E2

E1: bandas personales y barrido de la tolerancia. E2: competidores que se abstienen y Tabla 3 con ± DE (sección al final). E3 (reproducibilidad) no se ha hecho todavía.

## E1: bandas personales y barrido de la tolerancia

## Diseño, con las dos decisiones que se apartan de la especificación

- **Cohorte.** `cohort_labs12_v2.parquet` (paso 3 de `build_mimic_cohort.py`). Tiene las mismas 55 600 admisiones y añade, por examen:
  - la dispersión personal d = 1,4826 × MAD de los valores previos (los mismos que dan la línea base);
  - la re-medición r: el primer valor estrictamente posterior al del panel y a no más de 24 h.

  Los valores del panel y las líneas base recalculados coinciden exactamente con los del paso 2.
- **Banda.** La banda es [C − d, C + d], con C = b + resolución/4, el mismo umbral de `run_mimic.py`. Por tanto s = (V − C)/d.
  - La especificación decía centrar en b. Centrar en C garantiza que S_t y Neg_t quedan del lado correcto del umbral. Así la verdad observada sigue siendo una completación para todo t, y la diferencia es de un cuarto de la resolución.
  - Piso de d: un cuarto de la resolución, que afecta al 3,8 % de las celdas. Por examen va de 0,2 % a 0,5 % en plaquetas, hematocrito y leucocitos, y llega a 20,8 % en creatinina, por redondeo de los valores previos.
  - Variante `per_test`: un δ_m por examen, la mediana de d sobre la cohorte.
- **Muestra.** 30 726 admisiones (55 %) tienen re-medición de los doce exámenes. Las corridas muestrean solo de ese subconjunto, así que la verdad por re-medición es un contexto completo y no hay celdas híbridas. En consecuencia, la población difiere de la de la Tabla 3.
- **Parámetros.** N = 1 000, 10 repeticiones. Configuraciones: MCAR 10 %, MNARc 10 % y la referencia 0 %. t ∈ {0,25; 0,5; 1; 1,5; 2; 3}. 30 corridas, 105 s.
- **Verificación.**
  - `test_rho_int.py` compara ρ_int con la fuerza bruta sobre todas las completaciones (9 038 casos) y con `ptb.rho_int`.
  - `verify_regression.py`: ALL CHECKS PASSED.
  - Bajo la verdad observada, todas las corridas cumplen las aserciones de garantía.

## Números en el orden en que irían en el paper

### §7.2, datos clínicos: la re-medición como verdad empírica (`mimic_remeasure_stats.csv`)

- De las celdas del panel, 74,8 % tienen re-medición dentro de 24 h, con mediana de 12,7 h.
- Mediana de |s − s_re|: 0,67 semianchos (bandas personales). Cuantil 90: 3,4. Cuantil 99: 16.
- Fracción de celdas decididas (|s| ≥ t) que la re-medición pasa al otro lado de C, con bandas personales:

  | t | 0,25 | 0,5 | 1 | 1,5 | 2 | 3 |
  |---|---|---|---|---|---|---|
  | cambio de lado | 20,5 % | 18,4 % | 15,1 % | 13,2 % | 12,0 % | 10,8 % |

  Con una banda por examen, la tasa en t = 3 es de 7,0 %.
- **Esta tasa no cae a cero.** En un análisis exploratorio (no guardado en CSV), en t = 10 sigue entre 5 % y 12 %, según el tiempo entre mediciones. Varía mucho por examen en t = 3: potasio 25 %, glucosa 23 %, BUN 4 %, plaquetas 4 %.
- La lectura es que una segunda medición dentro de 24 h registra sobre todo cambio clínico (tratamiento) y repeticiones pedidas porque el primer valor era anómalo, no error de medición. Por eso la hipótesis de la Prop. 7.1 ("exacto dentro de la tolerancia") **no se cumple con la re-medición para ningún t del barrido**. El paper debe decirlo así. La especificación suponía lo contrario.

### §7.3, barrido de la tolerancia (`mimic_T_sweep_queries.csv`, figura `fig_tsweep`)

Preguntas paciente–patrón (patrones verdaderos de hasta tres exámenes con soporte ≥ 5 %), MCAR 10 %, bandas personales. Media de 10 repeticiones; la DE de cada cifra es ≤ 0,01.

| t | abstención | error de la certificación donde decide (vs. re-medición) | MICE en las preguntas abiertas (vs. re-medición) |
|---|---|---|---|
| 0,25 | 17,6 % | 10,5 % | 36,9 % |
| 0,5 | 26,5 % | 8,6 % | 33,5 % |
| 1 | 44,7 % | 5,9 % | 26,7 % |
| 1,5 | 58,9 % | 4,6 % | 22,6 % |
| 2 | 69,3 % | 3,7 % | 20,3 % |
| 3 | 81,7 % | 3,0 % | 17,9 % |

- Con una banda por examen, el error de la certificación donde decide en t = 3 es de 1,9 %, con abstención de 83,6 %.
- Con MNARc al 10 %, las cifras son casi iguales: 3,0 % de error y 80,8 % de abstención en t = 3.
- **Contra la verdad observada**, la certificación no yerra en ninguna corrida, ni con bandas. La abstención va de 17,6 % (t = 0,25) a 81,7 % (t = 3).
- Error de los métodos puntuales sobre todas las preguntas, MCAR 10 %:

  | método | vs. observada | vs. re-medición |
  |---|---|---|
  | MICE | 2,3 % | 15,2 % |
  | kNN | 3,0 % | 15,2 % |
  | media | 3,3 % | 15,5 % |
  | mediana | 3,4 % | 15,2 % |
  | ausente | 4,1 % | 13,8 % |

  Sin faltantes (0 %), el panel completo yerra 14,3 % contra la re-medición. Esa cifra es el piso que impone la deriva clínica.

### §7.3, ρ_int como predictor (`mimic_T_rho_int.csv`, panel derecho de `fig_tsweep`)

Intensiones certificadas en t = 0,25, frecuentes en el contexto superior, agrupadas en deciles de ρ_int. La tabla da la fracción que son intensiones del contexto re-medido (MCAR 10 %, bandas personales):

| decil | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| ρ_int medio | 0,26 | 0,30 | 0,36 | 0,43 | 0,51 | 0,57 | 0,70 | 0,88 | 1,16 | 2,10 |
| precisión | 91,5 % | 94,3 % | 96,5 % | 98,0 % | 97,1 % | 99,2 % | 99,6 % | 99,7 % | 99,9 % | 100 % |

- El patrón se repite con MNARc, sin faltantes y con una banda por examen: de 89–94 % en el primer decil a 100 % en el décimo.
- Precisión de las intensiones certificadas contra la re-medición: 97,6 % en t = 0,25, 99,9 % en t = 1 y 100 % desde t = 1,5.
- Recall contra la re-medición: 90,9 % en t = 0,25, 45,2 % en t = 1 y 3,3 % en t = 3.

### Conceptos exactos

- La certificación no recupera ningún concepto exacto en ningún t, banda o verdad (`rec_con` = 0). Las bandas no cambian lo que el paper ya dice al 10 % MCAR con δ = 0.

### Resumen (abstract)

- Hoy dice "certification never erred, as the theory guarantees for exact observed values". Sigue siendo cierto contra la verdad observada.
- Lo nuevo y defendible: "against a second measurement within 24 hours, the error of certification fell from 10.5% to 3.0% as t rose from 0.25 to 3, while abstention rose from 18% to 82%; the robustness ρ_int ranked certified intents from 91.5% to 100% correct."
- Hay que elegir qué cifras salen, porque tu guía permite tres o cuatro en el resumen.

## Frases del manuscrito que habría que cambiar

Las líneas se refieren a `revision_nov2026/personalized_thresholds_IJAR.tex` tal como está hoy.

| Línea | Qué dice hoy | Qué cambia |
|---|---|---|
| 71 | Resumen: "certification never erred, as the theory guarantees for exact observed values; with 10% … left 8% … 28%" | Añadir el resultado del barrido y de ρ_int; decidir qué cifras se quedan |
| 105 | Viñeta *Missing values* de la introducción: "certification never erred, as Proposition 7.1 guarantees" | Añadir que contra la re-medición el error baja con t y que ρ_int lo ordena |
| 706 | *Clinical data*: describe solo la mediana personal | Añadir la dispersión personal d (MAD), las bandas y la re-medición en 24 h, con el subconjunto de 30 726 admisiones |
| 712 | *Measurement noise*: escenario B con ruido < δ_m | Relegar o retirar, diciendo que se cumple por construcción; el barrido lo sustituye |
| 714 | *What is measured* | Definir el error de la certificación donde decide y la precisión por decil de ρ_int |
| 718 | *Guarantees*: "121 distinct settings … 780 … 430" | Sumar las 30 corridas del barrido, o dejarlas aparte para no mezclar recuentos |
| 753 | *Clinical data, patients* | Párrafo nuevo *Tolerance sweep* después de este, con `fig_tsweep` y las tablas de arriba |
| 762 | *Measurement noise* (resultados) | Coherente con la decisión de la l. 712 |
| 819 | Conclusión, "Two practical rules follow…" | La regla "report the concepts of large ρ" tiene ahora apoyo empírico (deciles de ρ_int) |
| nuevo §7.4 | — | Párrafo que distinga lo medido con δ = 0 (certificación de Burmeister–Holzer) de lo medido con bandas, y que diga que la re-medición en 24 h no satisface la hipótesis de la Prop. 7.1 |


---

## E2: competidores que se abstienen (`run_abstain.py`, `mimic_abstain.csv`, figura `fig_risk_coverage`)

### Diseño

- **Competidores.** Hay dos, y cada uno usa 20 tablas completadas:
  - imputación múltiple con `IterativeImputer(sample_posterior=True)`, con las líneas base como predictores;
  - kNN (k = 5) con donantes tomados de un remuestreo bootstrap.
- **Voto.** Para cada pregunta paciente–patrón, el voto es la fracción de tablas que responden «sí». Al nivel α, el competidor responde «sí» si el voto es ≥ 1 − α, «no» si es ≤ α, y se abstiene en medio. Se probaron α ∈ {0; 0,05; 0,1; 0,2; 0,3; 0,5}. α = 0 exige unanimidad, y α = 0,5 es la mayoría, sin abstención.
- **Diseño «main».** Reutiliza las muestras y máscaras de `run_mimic.py` con N = 1 000 y la verdad es el panel completo. Reproduce exactamente los números publicados: 8,0 % de preguntas abiertas, y MICE con 2,26 % de error en todas y 28,1 % en las abiertas.
- **Diseño «sweep».** Usa las muestras de E1, con las verdades obs y re-medición.
- **Configuraciones.** MCAR y MNARc al 10 %, 10 repeticiones, 40 corridas, unos 15 min.

### Resultados (MCAR 10 %; media de 10 repeticiones, DE ≤ 0,4 puntos)

**Verdad = panel completo (diseño main).** La certificación responde 92,0 % de las preguntas sin ningún error.

| competidor | α | cobertura | riesgo (error entre las respondidas) | responde, de las abiertas por la certificación | error en esas |
|---|---|---|---|---|---|
| imputación múltiple | 0 (unanimidad) | 93,3 % | 0,04 % | 15,9 % | 3,0 % |
| imputación múltiple | 0,1 | 94,1 % | 0,15 % | 26,7 % | 6,4 % |
| imputación múltiple | 0,3 | 96,7 % | 0,88 % | 58,4 % | 18,3 % |
| imputación múltiple | 0,5 (mayoría) | 100 % | 2,31 % | 100 % | 28,8 % |
| kNN bootstrap | 0 | 94,8 % | 0,74 % | 34,9 % | 25,1 % |
| MICE (una imputación) | — | 100 % | 2,26 % | 100 % | 28,1 % |
| ausente | — | 100 % | 3,81 % | 100 % | 47,4 % |

- **La imputación múltiple con unanimidad es un competidor serio.** Responde algo más que la certificación (93,3 % frente a 92,0 %) con un error casi nulo (0,04 %). Aun así, ese error no es cero, y no hay ninguna garantía que lo acote.
- Con MNARc los resultados son parecidos:
  - imputación múltiple con α = 0: 92,8 % de cobertura y 0,05 % de riesgo;
  - en las abiertas: 6,2 % de error.
- **La línea base «ausente»** es la peor en las preguntas abiertas: yerra en 47,4 % de ellas, casi al azar.

**Verdad = re-medición (diseño sweep).** La imputación múltiple no puede bajar de 93,0 % de cobertura, porque solo se abstiene cuando sus tablas discrepan, y su riesgo allí es 12,8 %. La certificación con bandas recorre toda la curva:

| t | cobertura | riesgo |
|---|---|---|
| 0,25 | 82,4 % | 10,5 % |
| 1 | 55,3 % | 5,9 % |
| 3 | 18,3 % | 3,0 % |

Con una banda por examen, en t = 3 la cobertura es 16,4 % y el riesgo 1,9 %. Ningún competidor ofrece un punto de operación con riesgo bajo frente a la re-medición.

### Tabla 3 regenerada (`make_table_results.py` → `table3.tex`, `table3_datasets.tex`)

- Media ± DE sobre las 10 repeticiones. Todos los valores coinciden con `c2_tabla3_propuesta.md`.
- Fila nueva «Absent (least completion)», con las celdas de preguntas que antes faltaban:
  - 3,8 ± 0,1 % de error en todas las preguntas;
  - 47,4 ± 1,1 % en las abiertas.
- Hay una tabla suplementaria por dataset público. La precisión de la media va de 0,40 (breast cancer) a 0,96 (diabetes).
- `summarize_results.py` cuenta ahora las configuraciones distintas: 121 (78 públicas y 43 clínicas), con 780 y 430 corridas, como dice el paper.

### Qué cambia en el paper si se integra E2

- **§7.2, *Methods*:** describir los dos competidores que se abstienen.
- **§7.3:**
  - un párrafo *Abstaining competitors* con la figura `fig_risk_coverage`;
  - la Tabla 3 reemplazada por `table3.tex`, y `table3_datasets.tex` como tabla suplementaria.
- **§7.4:** la ventaja de la certificación frente a la re-medición está en que gradúa la abstención y la acota con una garantía, no en el número de preguntas que responde. Con la verdad completa, la imputación múltiple unánime responde un poco más, casi sin errores, pero sin garantía.
