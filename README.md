# Fluorescence microarray 분석 — Phase 1

## 더블클릭으로 실행

`분석 실행.bat`를 더블클릭하면 `input` 폴더에서 TIFF 파일을 선택하는 창이 열립니다.
파일을 고른 뒤 **예**를 누르면 자동 검출, **아니요**를 누르면 직접 네 모서리를 선택합니다.
수동 선택은 왼쪽 위 → 오른쪽 위 → 왼쪽 아래 → 오른쪽 아래 순서이며, 확인 후 **Analyze**를 누릅니다.
결과는 `results/입력파일명_실행시각/`에 저장되고, 완료 창에 경로가 표시됩니다.
Python과 `requirements.txt`의 패키지는 먼저 설치되어 있어야 합니다.

## Phase 2 자동 grid 검출

Phase 1의 raw 로딩, 11×11 ROI와 통계/배경/CSV 측정 함수는 그대로 사용합니다. `--grid`를 생략하면 기존 수동 모드로 실행됩니다.

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe main.py chip01.tif --grid auto --output results
.\.venv\Scripts\python.exe main.py chip01.tif --grid manual --output results
```

자동 모드는 다음 순서로 동작합니다.

1. 검출 전용 float 배열에 log 변환과 두 Gaussian scale의 차이를 적용해 완만한 배경과 밝기 차이를 줄입니다. 이 배열은 intensity 계산에 쓰지 않습니다. 국소 최대값을 후보로 추출하고 작은 주변 영역의 중심을 계산합니다.
2. 최근접 후보 사이 displacement의 2D histogram에서 반복되는 두 방향의 pitch/rotation/skew를 추정합니다. 모든 후보가 동일하게 한 표씩 기여하므로 밝은 먼지 하나가 grid를 결정하지 않습니다.
3. 격자 좌표계에서 정수 간격에 가까운 이웃을 연결해 block 후보를 분리합니다. 각 block의 pitch를 다시 추정하고 phase와 가장 많이 지지되는 35×35 index window를 찾습니다.
4. bilinear 전역 geometry를 `soft_l1` robust least squares로 fitting하고, pitch의 일정 비율 안에 있는 후보만 재매칭하여 outlier를 제외합니다. 각 grid node는 개별 밝기 최대 위치로 이동하지 않습니다.
5. 충분한 행·열 지지가 있으면 전체 grid가 공유하는 낮은 차수의 cubic residual field를 fitting합니다. 최대 보정이 pitch의 15% 이내일 때만 적용합니다. 개별 ROI에 대한 intensity 최적화는 하지 않습니다.
6. 35개 행·열 각각의 후보 지지도, 후보 매칭률, 잔차, ROI 경계/겹침으로 QC 상태를 결정합니다. `OK` grid에 기존 Phase 1 측정을 적용합니다.

Robust fitting 구현은 [SciPy least_squares 공식 문서](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)의 `soft_l1` loss를 사용합니다. Threshold는 후보 생성에만 사용하며 최종 grid는 밝기순으로 결정하지 않습니다. 기본 confidence 기준은 매칭률 ≥60%, 각 행·열에서 ≥3개 후보, RMS < 최소 pitch의 12%입니다. 이 수치는 검출 지지도 기준이며 확률 또는 accuracy 보증은 아닙니다.

현재 검출은 x/y 축에서 약 ±25° 이내의 회전, 약한 skew, 11픽셀보다 큰 pitch를 대상으로 합니다. 완전히 신호가 없는 가장자리 행/열에서는 35×35 시작 위치가 모호해질 수 있으며, 심한 변형·넓은 오염·여러 pitch가 섞인 영상에서는 자동 검출을 신뢰하기 어렵습니다. 자동 결과 QC를 반드시 확인하십시오. 보정은 geometry 공유로 편향을 줄이지만, 후보 검출 자체가 signal에 의존하므로 편향이 완전히 없다고 보장하지 않습니다.

### 여러 block과 출력

한 TIFF에 여러 block이 있으면 검출된 block을 위→아래, 같은 높이에서는 왼쪽→오른쪽으로 번호 매깁니다. **번호는 검출된 후보의 순서**이며, 검출되지 않은 block을 포함한 실험 layout ID가 아닙니다.

```powershell
.\.venv\Scripts\python.exe main.py "input/YE chip1.tif" --grid auto --block 7 --output results
```

자동 출력은 기존 수동 결과를 보호하도록 `results/auto/` 아래에 저장합니다.

```text
results/auto/
    chip01_auto_grid_QC.png       # 모든 block 후보: 초록 OK, 주황 낮은 confidence/오류
    chip01_detection.json         # block별 상태, pitch, 회전, skew, center, 잔차
    block01/
        chip01_auto_grid_QC.png   # 보정된 자동 ROI 후보
        chip01_grid.json          # 자동 grid의 실수 좌표와 검출 정보
        chip01_confidence.csv     # spot별 후보 거리와 QC, 측정 CSV와 분리
        chip01_results.csv        # 기존 Phase 1과 동일한 열
        chip01_QC.png             # 실제 측정에 사용한 동일한 ROI bounds
        chip01_measurement.json   # auto/manual_fallback 출처와 실제 grid
```

자동 grid JSON은 실수 중심을 저장하고, 측정 CSV와 QC는 Phase 1과 동일한 반올림과 정확히 121픽셀 ROI를 사용합니다. 자동 후보 QC와 실제 측정 QC는 수동 fallback을 사용한 경우 서로 다를 수 있습니다. 실제 측정에는 항상 `*_QC.png`가 대응합니다.

같은 출력 경로를 재사용하면 같은 이름의 생성 파일을 덮어씁니다. 검출 결과가 달라져 이번 실행에서 측정하지 않은 block의 이전 CSV가 남을 수 있으므로, **검증 run별로 다른 `--output` 폴더를 사용**하고 해당 실행의 detection/measurement JSON을 확인하십시오.

### 실패와 manual fallback

`LOW_CONFIDENCE`, `OUT_OF_BOUNDS`, `GRID_MISMATCH`이면 먼저 자동 QC와 상태를 저장하고, 기본값인 `--fallback manual`에서 Phase 1의 4점 선택 창을 엽니다. 해당 block을 선택해 분석하면 출처가 `manual_fallback`으로 기록되고 클릭 좌표도 저장됩니다. 창을 닫으면 취소합니다.

GUI 없이 자동 결과만 검사하려면:

```powershell
.\.venv\Scripts\python.exe main.py chip01.tif --grid auto --fallback stop --output results/check01
```

이 경우 낮은 confidence block은 QC만 저장하고 intensity 결과를 생성하지 않습니다. 다른 `OK` block은 측정하며, 누락/실패가 있음을 알리기 위해 종료코드 1을 반환합니다. 후보 자체가 없으면 상태 JSON과 표시 영상만 저장합니다. 이후 `--grid manual`로 실행할 수 있습니다.

### Manual vs automatic 비교

기존 Phase 1에서 저장한 `*_corners.json`을 지정합니다. 같은 raw TIFF에서 두 grid를 **동일한 측정 함수**로 계산합니다. 여러 block이 있으면 manual grid 중심에 가장 가까운 block과 비교합니다. reference가 다른 block에 대응하지 않는지 QC에서 확인하십시오.

```powershell
.\.venv\Scripts\python.exe main.py "input/YE chip1.tif" --grid auto --fallback stop --compare-manual "results/YE chip1_corners.json" --output results/compare01
```

대응 block 폴더에 `*_comparison.csv`(실수 중심 차이, grid displacement, 실제 정수 ROI 중심 차이, Mean/Median 및 차이)와 `*_comparison_summary.json`(평균·중앙·최대 변위와 intensity 절대 차이)을 저장합니다. 낮은 confidence grid도 경계·겹침 조건이 유효하면 진단용 비교를 저장할 수 있지만, 이것은 자동 측정 승인을 의미하지 않습니다. 좌표 기준은 모두 0-based입니다.

### 테스트

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

Phase 1 회귀 검증에 더해 누락 spot·불균일 배경·밝은 artifact·회전/skew·매끄러운 변형의 합성 영상, blank 이미지 실패, 같은 grid 비교값이 0인지 확인합니다. 실제 이미지에 대한 ImageJ 일치 여부는 별도의 실험 검증이 필요합니다.

실제 `input/YE chip1.tif`의 초기 Phase 2 검증 결과는 `results/phase2_validation/auto/`에 저장했습니다. 16개 block 후보 중 10개는 OK, 6개(2, 7, 9, 10, 13, 15)는 LOW_CONFIDENCE입니다. 기존 수동 좌표에 대응하는 7번 block의 평균/최대 grid displacement는 2.630/5.004픽셀입니다. 해당 block의 자동 intensity CSV는 생성하지 않았고 진단용 비교만 저장했습니다. 이는 완전 자동 검출이 모든 block에서 확정된 상태가 아님을 의미합니다.

사용자가 선택한 TL, TR, BL, BR spot 중심을 bilinear interpolation하여 35×35 grid를 생성하고, 각 spot의 **11×11 raw TIFF 픽셀**을 측정합니다. 자동 검출, grid correction, spot 이동은 하지 않습니다.

## 설치와 실행 (Windows PowerShell)

Python 3.9 이상과 데스크톱 GUI 환경이 필요합니다.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe main.py chip01.tif --output results
```

이미 설치된 Python 환경에서는 다음과 같이 실행할 수 있습니다.

```powershell
python -m pip install -r requirements.txt
python main.py chip01.tif --output results
```

1. **TL → TR → BL → BR** 순서로 네 corner spot의 중심을 왼쪽 클릭합니다. 영상의 x는 오른쪽, y는 아래쪽으로 증가합니다.
2. **마우스 휠**로 커서 위치를 중심으로 확대/축소하고, **오른쪽 버튼 드래그**로 이동합니다. **Fit image**는 전체 보기로 돌아갑니다. 확대/이동 중에도 선택한 좌표는 유지됩니다. 도구 모음의 확대/이동을 사용할 경우에는 클릭하기 전 해당 모드를 해제합니다.
3. 잘못 클릭했으면 **Undo**(마지막 점 삭제) 또는 **Reset**(전체 삭제)을 누릅니다.
4. 네 점을 선택하면 1225개 ROI가 초록색으로 표시됩니다. 확대하여 위치를 확인합니다.
5. **Analyze**를 누르면 분석하고 파일을 저장합니다. 창을 그냥 닫으면 취소됩니다.

입력은 단일 2D grayscale TIFF입니다. 여러 페이지/채널의 stack과 RGB는 임의로 변환하지 않고 오류로 처리합니다. TIFF dtype은 그대로 유지하며 raw 배열은 읽기 전용입니다. 표시용 8-bit 영상은 별도 배열에 1–99.8 percentile contrast stretch를 적용합니다. 표시용 영상은 측정에 사용하지 않습니다.

## 출력

```text
results/
    chip01_results.csv
    chip01_QC.png
    chip01_corners.json
```

CSV는 1225행과 아래 열을 갖습니다. UTF-8 BOM을 포함하여 Excel에서도 열 수 있습니다.

```text
Spot,Row,Column,Center_X,Center_Y,Mean,Median,SD,CV,Min,Max,Total_Intensity,
BG_Mean,BG_Median,BG_SD,BG_CV,Pct_BG_1SD,Pct_BG_2SD
```

JSON에는 클릭한 네 점의 실수 좌표, 입력 이미지 정보, ROI 및 배경 설정을 저장합니다. 같은 이름으로 재실행하면 기존 출력 파일을 덮어씁니다.

## 좌표와 통계 정의

- Spot ID는 행 우선으로 `R01C01`부터 `R35C35`까지 부여합니다. Row/Column은 1부터 시작합니다.
- 픽셀 중심 좌표는 **0부터 시작**합니다. 보간된 중심은 `floor(coordinate + 0.5)`로 반올림합니다. CSV의 Center_X/Y는 실제 분석에 사용한 정수 중심입니다. 원래 보간 좌표는 저장된 네 corner로 재생성할 수 있습니다.
- 중심 `(x,y)`의 ROI는 `raw[y-5:y+6, x-5:x+6]`이며 정확히 121픽셀입니다. 경계를 벗어나거나 다른 foreground ROI와 겹치면 분석하지 않고 오류를 표시합니다.
- Mean/Median/Min/Max/Total_Intensity는 원본 강도의 통계입니다. 계산 시 overflow 방지를 위해 float64로 연산하며 intensity normalization은 하지 않습니다. 매우 큰 64-bit 정수(2^53 초과)는 float64 정밀도 한계가 있으므로 일반적인 8/16/32-bit fluorescence TIFF를 권장합니다.
- **SD와 BG_SD는 표본 표준편차(`ddof=1`)**입니다. CV/BG_CV는 SD/Mean×100이며 평균이 0이면 `nan`을 기록합니다.
- threshold 비율은 foreground 원본값이 `BG_Mean + n×BG_SD`를 **엄격하게 초과(>)**하는 픽셀 수 / 121 × 100입니다. 배경 차감 intensity는 이번 출력에 포함하지 않습니다.

## Local background

기본 배경 후보는 실제 정수 중심을 기준으로 한 **21×21 정사각형**입니다. 이 영역에서 **1225개 모든 foreground ROI의 합집합을 제외**합니다. 따라서 자기 ROI뿐 아니라 이웃의 foreground ROI도 배경에 포함되지 않습니다. 영상 경계에서는 배경 영역만 잘라냅니다. 남은 배경 픽셀이 2개 미만이면 분석을 중단합니다.

배경 크기를 변경하려면 11보다 큰 홀수를 지정합니다.

```powershell
python main.py chip01.tif --output results --background-size 25
```

Spot의 형광 꼬리가 11×11 밖으로 퍼지면 해당 픽셀은 배경에 들어갈 수 있습니다. 실제 영상의 spot 크기와 배경 영역을 고려해 설정해야 합니다. 자동 배경 영역 보정은 하지 않습니다.

## QC와 ImageJ 검증

QC PNG는 원본 이미지와 동일한 가로·세로 픽셀 수로 저장합니다. 분석에 전달한 동일한 ROI bounds로 `(x-5,y-5)`부터 `(x+5,y+5)`까지 테두리를 그립니다. 화면 preview에서는 픽셀 경계에 맞춰 half-pixel 경계를 사용합니다. QC 영상은 표시 contrast를 적용한 영상이고 정량 데이터는 CSV에 있습니다.

ImageJ에서 **원본 TIFF**를 열고 CSV의 정수 중심 `(x,y)`에 대해 왼쪽 위 `(x-5,y-5)`, 폭/높이 `11`인 rectangle을 만듭니다. 측정 영역, intensity calibration, threshold 제한 등 ImageJ 측정 설정이 같은지 확인한 뒤 Mean, Median, StdDev, Min, Max, RawIntDen을 비교합니다. 이 프로그램은 intensity calibration을 적용하지 않습니다. 배경은 위의 exclusion mask 방식과 동일하게 맞춰 비교해야 합니다.

실제 실험 이미지와 기존 ImageJ 결과는 제공되지 않았으므로 그 데이터의 일치 여부는 별도로 확인해야 합니다.

## 함수와 검증

`main.py`에는 로딩, 표시 변환, 수동 선택, grid 생성, ROI 추출, foreground/background 통계, threshold 계산, CSV 및 QC 저장 함수가 분리되어 있습니다.

```powershell
python -m unittest discover -s tests -v
```

합성 uint16 TIFF로 raw 값 보존, 회전·skew 보간, 121픽셀 통계, 배경에서 이웃 ROI 제외, CSV 행 수 및 QC 좌표를 검증합니다. GUI 클릭은 실제 데스크톱에서 확인해야 합니다.
