# Network Time Series Forecasting using Spectral Graph Wavelet Transform

This repo contains the source code for network time series forecasting using SGWT. 

For the detail, please refer to our paper: [https://doi.org/10.1016/j.ijforecast.2023.08.006](https://doi.org/10.1016/j.ijforecast.2023.08.006).

## Description

- Code
  - `method.R` is a code for functions used for forecasting.
  - `seoulmetro.R` and `covidseoul.R` are codes for real data analysis.

- Data
  - `seoulmetro` contains data of daily number of people getting on and off the Seoul Metropolitan Subway in South Korea.
  - `covidseoul` contains data of daily number of newly confirmed cases in the Seoul Metropolitan Area.
  

## Code overview
We apply various network forecasting methods to analyze real-world data, Seoul Metropolitan Subway data and Seoul Metropolitan COVID-19 data.
The methods used in the analysis are

  - Nodewise ARIMA
  - LOCAAT-based method
  - GNARI
  - GNARI + LOCAAT
  - GFT-based method
  - SGWT-based method