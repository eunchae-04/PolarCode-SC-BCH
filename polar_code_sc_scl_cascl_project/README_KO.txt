Polar Code SC / SCL / CRC-aided SCL 기준선 프로젝트
====================================================

이 폴더의 polar_code_sc_scl_cascl_compare.c는 사용자가 첨부한 전체 C 소스입니다.
이 파일 안에 main(), Polar encoder, GA frozen-set 생성, AWGN 채널,
SC, SCL(L=2/4/8), CRC-aided SCL(L=2/4/8)이 이미 들어 있습니다.
SC/SCL/CA-SCL 디코더를 별도 파일로 새로 만들 필요가 없습니다.

중요:
- 현재 이 프로그램에는 BCH encoder/decoder가 들어 있지 않습니다.
- graph_guided_scflip.c/.h는 일반 SC 디코더와 아직 연결되지 않은 adapter 모듈이며,
  이 기준선 프로그램의 빌드 명령에 포함하지 않습니다.
- 따라서 현재 생성되는 그래프는 SC/SCL/CA-SCL 비교이며 BCH+SC 결과가 아닙니다.
- CA-SCL은 CRC-8(poly 0x07)을 사용합니다. BCH와 CRC는 다른 코드입니다.

[1] MSYS2 UCRT64에서 컴파일

프로젝트 폴더로 이동한 다음 실행:

gcc -std=c11 -O2 -Wall -Wextra -Wpedantic -o compare.exe polar_code_sc_scl_cascl_compare.c -lm

또는:

bash build_ucrt64.sh

[2] 짧은 동작 확인

./compare.exe -1 0 1 100 64 0.5 result_smoke

인수 의미:
  -1       시작 Eb/N0 (dB)
   0       종료 Eb/N0 (dB)
   1       Eb/N0 간격 (dB)
 100       각 디코더/SNR에서 관측할 비트 오류 목표 수
  64       Polar block length N
 0.5       Polar rate R (K=round(N*R))
 result_smoke 결과 저장 폴더

명령이 끝나면 result_smoke 폴더에 다음 파일이 생성됩니다.
  result_SC.txt
  result_SCL_L2.txt
  result_SCL_L4.txt
  result_SCL_L8.txt
  result_CASCL_L2.txt
  result_CASCL_L4.txt
  result_CASCL_L8.txt

각 결과 파일 열:
  EbNo_dB BER FER frames

[3] 연구 조건 실행 예시

N=256, R=1/3:
./compare.exe -1 3 0.5 1000 256 0.333333 result_N256_R033

N=512, R=1/2:
./compare.exe -1 3 0.5 1000 512 0.5 result_N512_R05

주의: 소스는 프레임 수 상한 없이 각 지점에서 목표 비트 오류 수에 도달할 때까지 실행합니다.
높은 Eb/N0에서 오류가 매우 드물면 한 지점이 오래 걸릴 수 있습니다.
먼저 목표 오류 수 50~100으로 동작 확인을 하고, 연구용 데이터는 충분한 오류 표본을
확보하도록 목표 수를 높이세요. 높은 SNR에서 오류가 전혀 발생하지 않으면 실행이
계속될 수 있으므로 진행 상황을 확인하고 필요하면 Ctrl+C로 중단하세요.

[4] MATLAB 그래프

MATLAB Current Folder를 이 폴더로 지정한 다음:

plot_compare_results('result_smoke')

연구 실행 결과라면 해당 결과 폴더 이름을 인수로 넣으세요.

[5] 현재 코드의 조건 및 해석 주의사항

- 지원 N: 64, 128, 256, 512, 1024 (2의 거듭제곱)
- K=round(N*R), K는 CRC_LEN=8보다 커야 합니다.
- Frozen set은 GA 방식으로 각 Eb/N0 지점에서 다시 계산됩니다.
- SC/SCL은 K개의 정보를 전송하고, CA-SCL은 K-8개의 payload와 CRC-8을 전송합니다.
  따라서 payload 유효 부호율이 서로 다르며, 결과를 비교할 때 이 오버헤드를 명시해야 합니다.
- 현재 출력 파일에는 BER, FER, 프레임 수만 기록됩니다. 평균 복호 시간은 아직 측정하지 않습니다.

BCH+SC 및 그래프 기반 후보 flip을 추가할 때는 BCH 코드 파라미터(길이, 정보 길이,
생성다항식 또는 정정 능력)와 교수님 그래프의 인덱스 기준(절대 u-index인지)을 확인한 뒤
별도의 알고리즘으로 추가해야 합니다. 임의의 BCH 파라미터를 적용한 결과를 교수님 방식의
성능 결과라고 해석하면 안 됩니다.
