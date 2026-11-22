# QUPLAS: proposte per il fascio di positronio

Pagina interattiva con quattro proposte per la linea di fascio di QUPLAS (positroni → Ps⁻ → fotodistacco laser → Ps). Si basa su un modello semplificato con i parametri di:

> A. Galanti et al., *Design of an experimental setup for continuous positronium beam formation via Ps⁻ production and laser photo-detachment*, Radiation Physics and Chemistry 251 (2027) 114388. doi:10.1016/j.radphyschem.2026.114388

## Contenuto

- `index.html`: la pagina, autonoma, da aprire nel browser o pubblicare con GitHub Pages.
- `quplas_bilancio_flusso.py`: lo stesso bilancio di flusso in Python (`python3 quplas_bilancio_flusso.py`).

## Le quattro proposte

1. **Misura finale su una sorgente ad alto flusso** (ad esempio NEPOMUC, continua e compatibile con la cavità CW).
2. **Contare gli elettroni del fotodistacco:** monitor in tempo reale della produzione di Ps, anello chiuso per la rigenerazione del convertitore e segnale di start per le coincidenze con il rivelatore finale.
3. **Convertitore in trasmissione, laser a pochi millimetri:** il Ps⁻ vive 479 ps, ogni millimetro conta.
4. **Specchi isolati e monitorati:** preparazione del convertitore in una camera separata, misura periodica dello spostamento dei modi trasversi.

## Limiti

I valori marcati "ipotesi" non vengono dall'articolo e vanno sostituiti con quelli reali. Il modello serve a confrontare scenari e non sostituisce le simulazioni complete della collaborazione.
