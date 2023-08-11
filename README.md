# BRCA_2023
P-DIAMOND BRAC Study

# 정량 단계
1. Run SiteQuantification.java
2. Run GISLog2Norm.java

# 정량에 사용하는 파일
- PSMs.txt
- QuanSpectra.txt

# 파일 별 처리 및 연결 관계

## 파일 별 관계

- QuanSpectra — <First Scan + Spectrum File> —> PSMs.txt

## 제외되는 경우

- QuanSpectra
    - NoQuanValues
    - RejectedByMethod
- PSMs.txt
    - Master proteins ≥ 2

## 파일 별 특징

### QuanSpectra

- Quan Info에서 정량 가능성이 기입되어 있음
    - 빈 칸: 정량 가능
    - NoQuanValues: 정량 불가능
    - RejectedByMethod: 정량 불가능
- PSM파일과 동일하게 채널 간 normalization이 적용되어 있음.
- 사이즈가 크기 때문에, 해당 파일에서는 정량 가능한 스펙트럼만 체크하는 역할로 사용

### PSMs.txt

- Quan Info에 정량 가능성이 기입되어 있어야 하지만.. 실제로는 빈 칸임.
따라서, QuanSpectra를 활용해서 얻어와야 함.
- Modifications에 modification site 정보가 있으나! PhosphoRS에서 correction된 값을 사용해야하기 때문에, “PhosphoRS: Best Site Probabilities” 칼럼을 사용해서 위치 정보를 확인해야 함.
- 펩타이드 서열과 정량 사이트를 활용해서 site 정보를 얻고, site 정보를 바탕으로 PeptideGroups.txt와 연결해서 사용.
- peptide sequence에 PTM에 따른 **대소문자 구분이 있음**.

### PeptideGroups.txt

- peptide sequence에 PTM에 따른 **대소문자 구분이 없음**.
- 또한, QuanInfo에서 NoQuanValues는 제외.
- 근데, PeptideGroup 정보를 나타내기 때문에, 사이트 정보가 뭉개지는 현상이 있음.

### ** PeptideGroups.txt에서 사이트 정보가 뭉개지는 현상 **

![PeptideGroups.txt](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/0c6cf106-6f23-4150-9ab1-bda5bebed5ec/Untitled.png)

PeptideGroups.txt

![PSMs.txt](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/221ba8ce-a275-4e65-bd1d-b3b55a7cbbd7/Untitled.png)

PSMs.txt

PeptideGroups.txt에서 phospho가 두 개인 경우를 보면 (두번째 줄), S5, S11로 되어 있음. 하지만 개별 PSM으로 보면, S5, S11이 두 개 있고, S5, S6가 하나 있음. 따라서, PeptideGroups.txt를 보는 것은 무의미할 듯.
