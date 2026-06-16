# AquaLNCRpred

A web server for predicting long non-coding RNA (lncRNA) sequences in *Eriocheir sinensis* (Chinese Mitten Crab) and related aquaculture organisms.

## 1. Download AquaLNCRpred

Download the source code from the GitHub repository:

[AquaLNCRpred GitHub Repository](https://github.com/malik010/AquaLNCRpred?utm_source=chatgpt.com)

### Option A: Download ZIP Archive

1. Open the repository link above.
2. Click **Code** → **Download ZIP**.
3. Extract the downloaded ZIP file to your desired directory.

### Option B: Clone Using Git

```bash
git clone https://github.com/malik010/AquaLNCRpred.git
cd AquaLNCRpred
```

## 2. Installation

Install the required R packages:

```r
install.packages(c(
  "shiny",
  "shinydashboard",
  "dashboardthemes",
  "markdown",
  "seqinr",
  "shinyjs",
  "protr",
  "randomForest",
  "caret",
  "tidyverse",
  "kernlab",
  "ftrCOOL",
  "doParallel",
  "xgboost"
))
```

## 3. Run AquaLNCRpred as a Local Web Server

After all required packages have been installed successfully:

1. Start R.

2. Load the `shiny` package:

```r
library(shiny)
```

3. Set the working directory to the AquaLNCRpred folder:

```r
setwd("path/to/AquaLNCRpred")
```

4. Run the application:

```r
runApp("app.R")
```

The AquaLNCRpred web interface will open automatically in your default web browser.
