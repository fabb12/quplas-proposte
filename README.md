# QUPLAS: proposte per il fascio di positronio

Pagina interattiva con tre proposte per la linea di fascio di QUPLAS (positroni → Ps⁻ → fotodistacco laser → Ps). Si basa su un modello semplificato con i parametri di:

> A. Galanti et al., *Design of an experimental setup for continuous positronium beam formation via Ps⁻ production and laser photo-detachment*, Radiation Physics and Chemistry 251 (2027) 114388. doi:10.1016/j.radphyschem.2026.114388

## Contenuto

- `index.html`: la pagina, autonoma, da aprire nel browser o pubblicare con GitHub Pages.
- `guida.html`: guida in nove sezioni (principio di equivalenza, difficoltà, catena e⁺→Ps⁻→Ps, interferometria, statistica, obiettivo, antigravità, interferometro Mach-Zehnder a laser del progetto reale, elenco delle fonti) con animazioni e calcolatori basati sugli articoli della collaborazione.
- `interferenza.html`: simulazione dell'interferenza di singoli positroni (QUPLAS-0, Como 2019), con l'animazione dell'intero sistema.
- `quplas_bilancio_flusso.py`: lo stesso bilancio di flusso in Python (`python3 quplas_bilancio_flusso.py`).

`index.html` include anche un'animazione schematica della linea di fascio (e⁺ → Ps⁻ → Ps → interferometro), con rese illustrative.

## Le tre proposte

1. **Fare la misura su una sorgente ad alto flusso** (ad esempio NEPOMUC, continua come la cavità laser).
2. **Convertitore in trasmissione, laser a pochi millimetri:** il Ps⁻ vive 479 ps, ogni millimetro conta.
3. **Proteggere e sorvegliare gli specchi:** preparazione del convertitore in una camera separata, misura periodica dello spostamento dei modi trasversi.

## Limiti

I valori marcati "ipotesi" non vengono dall'articolo e vanno sostituiti con quelli reali. Il modello serve a confrontare scenari e non sostituisce le simulazioni complete della collaborazione.
