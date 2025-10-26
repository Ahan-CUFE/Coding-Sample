This project explores whether investor sentiment extracted from Chinese social-media commentary can improve quantitative forecasts for the CSI 300 (HS300) index. A sequence of text-processing scripts generates sentiment features that are later fused with traditional market indicators to train stacked LSTM regressors. The workflow enables a direct comparison between purely financial baselines and sentiment-augmented models.

- Baseline Sentiment (1.1Only_snownlp.py):  Applies SnowNLP directly to each headline without any cleaning to produce onlySnownlp.csv.
- Cleaned Sentiment (1.2text_pre-processing.py): Converts Traditional Chinese to Simplified using OpenCC, removes non-text symbols, and rescored posts to create processed_text_and_snownlp.csv.
- POS-Aware Sentiment (1.3pre_jieba_snownlp.py): Adds Jieba segmentation with part-of-speech filtering before SnowNLP scoring, exporting pre_jieba_snownlp_avn_vn_saved.csv.
- Enhanced Sentiment (advantaged_nlp.py): Repeats cleaning and POS filtering on the Excel dataset with broader tag retention and anomaly handling, yielding advantaged_trained_nlp.csv.
- Timeseries Alignment (main.py): Aggregates multiple sentiment files to daily averages for visualization against market movements.

## Modeling Strategy
Three LSTM setups quantify the marginal value of sentiment:
 2.1lstm_traditional.py — Uses only historical market features (from hs300data.csv) with PCA compression and a 5-day lookback window.
 2.2noemotion_lstm.py — Mirrors the enhanced pipeline but intentionally omits the sentiment merge, isolating the effect of preprocessing.
 2.3advantaged_lstm.py — Merges daily sentiment (emotionvalue) with HS300 fundamentals, performs PCA to six components, and trains on 80% of observations with a 7-day window.

Each model scales inputs to [0, 1], reshapes them into sequences, and trains a two-layer LSTM with dropout regularization. Predictions are rescaled to original prices for evaluation and plotting.
