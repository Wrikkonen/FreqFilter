#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * ファイルの存在確認
 * ----
 * fname : ファイル名
 * f()   : 存在フラグを返す（存在するなら1、存在しないなら0）
 */
int isFile(char *fname)
{
	FILE *fp;
	
	fp = fopen(fname, "rb");
	if (fp) fclose(fp);
	return (fp)? 1: 0;
}

/*
 * 画像をファイルから読み込む
 * ----
 * fname         : ファイル名
 * rev, imv      : データバッファ（実部、虚部）
 * width, height : データサイズ（横幅、高さ）
 * scale, offset : 画像化情報（スケール、オフセット）
 * f()           : 成功フラグを返す（成功なら1、失敗なら0）
 */
int readRawFileEx(char *fname, double *rev, double *imv, int width, int height, int rWidth, int rHeight, int *scale, int *offset)
{
	int ret = 0, x, y, size, minv, maxv;
	unsigned short *imgv = NULL;
	FILE *fp;
	
	*scale = 0;
	for (;;) {
		// データ領域を確保する
		size = rWidth * rHeight;
		if ((imgv = (unsigned short *)malloc(sizeof(unsigned short)*size)) == NULL) break;
		// ファイルからピクセル値を読み込む
		if ((fp = fopen(fname, "rb")) == NULL) break;
		fread(imgv, sizeof(unsigned short), size, fp);
		fclose(fp);
		// 実部にピクセル値、虚部にゼロをセットする
		minv = maxv = (int)imgv[0];
		for (y = 0; y < rHeight; y++) {
			for (x = 0; x < rWidth; x++) {
				if (x < width && y < height) {
					if (imgv[width*y+x] < minv) minv = imgv[width*y+x];
					if (maxv < imgv[width*y+x]) maxv = imgv[width*y+x];
					rev[rWidth*y+x] = (double)imgv[width*y+x];
					imv[rWidth*y+x] = 0.0;
				} else {
					rev[rWidth*y+x] = 0.0;
					imv[rWidth*y+x] = 0.0;
				}
			}
		}
		if (maxv <= minv) break;
		*scale = maxv - minv;
		*offset = minv;
		// 正常終了
		ret = 1;
		break;
	}
	// データ領域を解放する
	if (imgv) free(imgv);
	return ret;
}

/*
 * 画像をファイルへ書き込む
 * ----
 * fname         : ファイル名
 * rev, imv      : データバッファ（実部、虚部）
 * width, height : データサイズ（横幅、高さ）
 * scale, offset : 画像化情報（スケール、オフセット）
 * f()           : 成功フラグを返す（成功なら1、失敗なら0）
 */
int writeRawFileEx(char *fname, double *rev, double *imv, int width, int height, int rWidth, int rHeight, int scale, int offset)
{
	int ret = 0, x, y, size;
	unsigned short *imgv = NULL;
	double maxv, minv;
	FILE *fp;
	
	for (;;) {
		// データ領域を確保する
		size = rWidth * rHeight;
		if ((imgv = (unsigned short *)malloc(sizeof(unsigned short)*size)) == NULL) break;
		// 実部にピクセル値、虚部にゼロをセットする
		maxv = minv = rev[0];
		for (y = 0; y < rHeight; y++) {
			for (x = 0; x < rWidth; x++) {
				if (x < width && y < height) {
					if (rev[rWidth*y+x] < minv) minv = rev[rWidth*y+x];
					if (maxv < rev[rWidth*y+x]) maxv = rev[rWidth*y+x];
				}
			}
		}
		if (maxv <= minv) break;
		// ピクセル値をスケールの範囲で正規化する
		for (y = 0; y < rHeight; y++) {
			for (x = 0; x < rWidth; x++) {
				if (x < width && y < height) {
					imgv[width*y+x] = (unsigned short)((rev[rWidth*y+x]-minv)/(maxv-minv)*scale+offset);
				}
			}
		}
		// ピクセル値をファイルに書き込む
		if ((fp = fopen(fname, "wb")) == NULL) break;
		fwrite(imgv, sizeof(unsigned short), size, fp);
		fclose(fp);
		// 正常終了
		ret = 1;
		break;
	}
	// データ領域を解放する
	if (imgv) free(imgv);
	return ret;
}

/*
 * 中間ファイル（パワースペクトル）を書き込む
 * ----
 * fname         : ファイル名
 * rev, imv      : データバッファ（実部、虚部）
 * width, height : データサイズ（横幅、高さ）
 * scale, offset : 画像化情報（スケール、オフセット）
 * f()           : 成功フラグを返す（成功なら1、失敗なら0）
 */
int writeRawFileSpect(char *fname, double *rev, double *imv,
					  int width, int height, int scale, int offset)
{
	int ret = 0, i, size, first = 1;
	unsigned short *imgv = NULL;
	double *psv = NULL, d, maxv, minv;
	FILE *fp;
	
	for (;;) {
		// データ領域を確保する
		size = width * height;
		if ((psv = (double *)malloc(sizeof(double)*size)) == NULL) break;
		if ((imgv = (unsigned short *)malloc(sizeof(unsigned short)*size)) == NULL) break;
		// パワースペクトラムを求める
		for (i = 0; i < size; i++) {
			d = sqrt(rev[i]*rev[i] + imv[i]*imv[i]);
			psv[i] = (d < 1.0E-4)? log(1.0E-4): log10(d);
			if (first) {maxv = minv = psv[i]; first = 0;}
			if (psv[i] < minv) minv = psv[i];
			if (maxv < psv[i]) maxv = psv[i];
		}
		// ピクセル値をスケールの範囲で正規化する
		for (i = 0; i < size; i++) {
			imgv[i] = (unsigned short)((psv[i]-minv)/(maxv-minv)*scale+offset);
		}
		// ピクセル値をファイルに書き込む
		if ((fp = fopen(fname, "wb")) == NULL) break;
		fwrite(imgv, sizeof(unsigned short), size, fp);
		fclose(fp);
		// 正常終了
		ret = 1;
		break;
	}
	// データ領域を解放する
	if (imgv) free(imgv);
	if (psv) free(psv);
	return ret;
}

/*
 * 変数の値を入れ替える
 * ----
 * a, b : 入れ替える変数
 */
void swap(double *a, double *b)
{
	double tmp = *a;
	*a = *b;
	*b = tmp;
}

/*
 * 低周波と高周波のデータを入れ替える
 * ----
 * imgv          : 画像バッファ
 * width, height : 画像サイズ（横幅、高さ）
 */
void replace2d(double *imgv, int width, int height)
{
	int x, y, xh, yh;
	
	xh = width / 2;
	yh = height / 2;
	for (y = 0; y < yh; y++) {
		for (x = 0; x < xh; x++) {
			swap(&imgv[width*y+x], &imgv[width*(y+yh)+x+xh]);
			swap(&imgv[width*(y+yh)+x], &imgv[width*y+x+xh]);
		}
	}
}

/*
 * 離散フーリエ変換を行う
 * ----
 * rev, imv : データバッファ（実部、虚部）
 * n        : データ数
 * inverse  : 逆変換フラグ（0:順変換、1:逆変換）
 * f()      : 成功フラグを返す（成功なら1、失敗なら0）
 */
int dft(double *rev, double *imv, int n, int inverse)
{
	int ret = 0, i, j;
	double *rew = NULL, *imw = NULL, theta, sign, d;
	
	for (;;) {
		// データ領域を確保する
		if ((rew = (double *)malloc(sizeof(double)*n)) == NULL) break;
		if ((imw = (double *)malloc(sizeof(double)*n)) == NULL) break;
		// 演算
		sign = (inverse)? 1.0: -1.0;
		for (j = 0; j < n; j++) {
			theta = 2.0 * M_PI * j / (double)n;
			rew[j] = imw[j] = 0.0;
			for (i = 0; i < n; i++) {
				rew[j] += rev[i] * cos(theta*i) - sign * imv[i] * sin(theta*i);
				imw[j] += imv[i] * cos(theta*i) + sign * rev[i] * sin(theta*i);
			}
		}
		// 結果を戻す（順変換の場合は1/nする）
		d = (inverse)? 1.0: (double)n;
		for (i = 0; i < n; i++) {
			rev[i] = rew[i] / d;
			imv[i] = imw[i] / d;
		}
		// 正常終了
		ret = 1;
		break;
	}
	// データ領域を解放する
	if (imw) free(imw);
	if (rew) free(rew);
	return ret;
}

/*
 * 高速フーリエ変換を行う
 * ----
 * rev, imv : データバッファ（実部、虚部）
 * n        : データ数
 * inverse  : 逆変換フラグ（0:順変換、1:逆変換）
 * f()      : 成功フラグを返す（成功なら1、失敗なら0）
 */
int fft(double *rev, double *imv, int n, int inverse)
{
	int i, j, k, m, n1, n2, t;
	double theta, wr, wi, dr, di;
	
	// バタフライ演算
	theta = ((inverse)? M_PI: -M_PI) / (double)n;
	for (i = n, m = 0; 2 <= i; i >>= 1, m++)
		;
	n2 = n;
	for (k = 0; k < m; k++) {
		theta *= 2.0;
		n1 = n2;
		n2 = n1 / 2;
		for (j = 0; j < n2; j++) {
			wr = cos(theta*j);
			wi = sin(theta*j);
			for (i = j; i < n; i+=n1) {
				t = i + n2;
				dr = rev[i] - rev[t];
				di = imv[i] - imv[t];
				rev[i] += rev[t];
				imv[i] += imv[t];
				rev[t] = dr*wr - di*wi;
				imv[t] = di*wr + dr*wi;
			}
		}
	}
	// データを入れ替える
	for (i = j = 0; i < n-1; i++, j+=n2) {
		if (i < j) {
			swap(&rev[i], &rev[j]);
			swap(&imv[i], &imv[j]);
		}
		for (n2 = n/2; n2 <= j; j-=n2, n2>>=1)
			;
	}
	// 順変換の場合は1/nする
	if (! inverse) {
		for (i = 0; i < n; i++) {
			rev[i] /= (double)n;
			imv[i] /= (double)n;
		}
	}
	return 1;
}

/*
 * 2次元のフーリエ変換を行う
 * ----
 * rev, imv      : データバッファ（実部、虚部）
 * width, height : 画像サイズ（横幅、高さ）
 * inverse       : 逆変換フラグ（0:順変換、1:逆変換）
 * f()           : 成功フラグを返す（成功なら1、失敗なら0）
 */
int ft2d(double *rev, double *imv, int width, int height, int inverse)
{
	int ret = 0, x, y, side;
	double *rew = NULL, *imw = NULL;
	
	for (;;) {
		// データの入れ替え
		if (inverse) {
			replace2d(rev, width, height);
			replace2d(imv, width, height);
		}
		// データ領域を確保する
		side = (width < height)? height: width;
		if ((rew = (double *)malloc(sizeof(double)*side)) == NULL) break;
		if ((imw = (double *)malloc(sizeof(double)*side)) == NULL) break;
		// X方向のフーリエ変換を行う
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				rew[x] = rev[width*y+x];
				imw[x] = imv[width*y+x];
			}
			fft(rew, imw, width, inverse);
			for (x = 0; x < width; x++) {
				rev[width*y+x] = rew[x];
				imv[width*y+x] = imw[x];
			}
		}
		// Y方向のフーリエ変換を行う
		for (x = 0; x < width; x++) {
			for (y = 0; y < height; y++) {
				rew[y] = rev[width*y+x];
				imw[y] = imv[width*y+x];
			}
			fft(rew, imw, height, inverse);
			for (y = 0; y < height; y++) {
				rev[width*y+x] = rew[y];
				imv[width*y+x] = imw[y];
			}
		}
		// データの入れ替え
		if (! inverse) {
			replace2d(rev, width, height);
			replace2d(imv, width, height);
		}
		// 正常終了
		ret = 1;
		break;
	}
	// データ領域を解放する
	if (imw) free(imw);
	if (rew) free(rew);
	return ret;
}

/*
 * フィルタ処理する
 * ----
 * rev, imv      : データバッファ（実部、虚部）
 * width, height : 画像サイズ（横幅、高さ）
 * low, high     : 周波数（下限、上限：0:未設定、単位：サイクル／ピクセル）
 * f()           : 成功フラグを返す（成功なら1、失敗なら0）
 */
int filterFreq(double *rev, double *imv, int width, int height, int low, int high)
{
	int x, y, r, xh, yh, index, rLow, rHigh;
	
	// 画像の中心に円をかいてフィルタ処理する
	xh = width / 2;
	yh = height / 2;
	rLow = low * low;
	rHigh = high * high;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			r = (x-xh)*(x-xh) + (y-yh)*(y-yh);
			index = width * y + x;
			// バンドパスフィルタ
			if (0 < low && 0 < high) {
				if (r < rLow || rHigh < r) rev[index] = imv[index] = 0.0;
			}
			// ハイパスフィルタ
			else if (low <= 0 && 0 < high) {
				if (r < rHigh) rev[index] = imv[index] = 0.0;
			}
			// ローパスフィルタ
			else if (0 < low && high <= 0) {
				if (rLow < r) rev[index] = imv[index] = 0.0;
			}
		}
	}
	return 1;
}

/*
 * 2のべき乗に切り上げる
 * ----
 * side: 辺の長さ (1以上)
 * f() : 結果
 */
int roundupBin(int side)
{
	int value, shift;
	
	if (side <= 0) return 0;
	shift = 0;
	for (value = side-1; value; value>>=1)
		shift++;
	return 1<<shift;
}

/*
 * 画像処理する
 * ----
 * dst, src      : ファイル名（出力、入力）
 * width, height : 画像サイズ（横幅、高さ）
 * f()           : 成功フラグを返す（成功なら1、失敗なら0）
 */
int processA(char *dst, char *src, int width, int height)
{
	char *tmp = "temp.raw";
	enum {forward, inverse};
	int ret = 0, rWidth, rHeight, size, scale, offset;
	double *rev = NULL, *imv = NULL;
	
	for (;;) {
		// データ領域を確保する
		rWidth = roundupBin(width);
		rHeight = roundupBin(height);
		size = rWidth * rHeight;
		if ((rev = (double *)malloc(sizeof(double)*size)) == NULL) break;
		if ((imv = (double *)malloc(sizeof(double)*size)) == NULL) break;
		// フィルタ処理
		if (readRawFileEx(src, rev, imv, width, height, rWidth, rHeight, &scale, &offset) == 0) break;
		if (ft2d(rev, imv, rWidth, rHeight, forward) == 0) break;
		if (writeRawFileSpect(tmp, rev, imv, rWidth, rHeight, scale, offset) == 0) break;
		if (ft2d(rev, imv, rWidth, rHeight, inverse) == 0) break;
		if (writeRawFileEx(dst, rev, imv, width, height, rWidth, rHeight, scale, offset) == 0) break;
		// 正常終了
		ret = 1;
		break;
	}
	// データ領域を解放する
	if (imv) free(imv);
	if (rev) free(rev);
	return ret;
}

/*
 * 画像処理する
 * ----
 * dst, src      : ファイル名（出力、入力）
 * width, height : 画像サイズ（横幅、高さ）
 * low, high     : 周波数（下限、上限：0:未設定、単位：サイクル／ピクセル）
 * f()           : 成功フラグを返す（成功なら1、失敗なら0）
 */
int processB(char *dst, char *src, int width, int height, int low, int high)
{
	char *tmp = "temp.raw";
	enum {forward, inverse};
	int ret = 0, rWidth, rHeight, size, scale, offset;
	double *rev = NULL, *imv = NULL;
	
	for (;;) {
		// データ領域を確保する
		rWidth = roundupBin(width);
		rHeight = roundupBin(height);
		size = rWidth * rHeight;
		if ((rev = (double *)malloc(sizeof(double)*size)) == NULL) break;
		if ((imv = (double *)malloc(sizeof(double)*size)) == NULL) break;
		// フィルタ処理
		if (readRawFileEx(src, rev, imv, width, height, rWidth, rHeight, &scale, &offset) == 0) break;
		if (ft2d(rev, imv, rWidth, rHeight, forward) == 0) break;
		if (filterFreq(rev, imv, rWidth, rHeight, low, high) == 0) break;
		if (writeRawFileSpect(tmp, rev, imv, rWidth, rHeight, scale, offset) == 0) break;
		if (ft2d(rev, imv, rWidth, rHeight, inverse) == 0) break;
		if (writeRawFileEx(dst, rev, imv, width, height, rWidth, rHeight, scale, offset) == 0) break;
		// 正常終了
		ret = 1;
		break;
	}
	// データ領域を解放する
	if (imv) free(imv);
	if (rev) free(rev);
	return ret;
}

int main(void)
{
	char dst[256], src[256];
	int flag, width, height, low, high;
	
	printf("入力ファイル： ");				 scanf("%s", src);
	if (isFile(src) == 0) {
		printf("ファイル[%s]が見つかりません。\n", src);
		return 0;
	}
	printf("横幅（ピクセル）： ");			 scanf("%d", &width);
	printf("高さ（ピクセル）： ");			 scanf("%d", &height);
	printf("低周波（サイクル／ピクセル）： ");	 scanf("%d", &low);
	printf("高周波（サイクル／ピクセル）： ");	 scanf("%d", &high);
	printf("出力ファイル： ");				 scanf("%s", dst);
	
	flag = processB(dst, src, width, height, low, high);
	printf((flag)? ">成功しました。\n": ">失敗しました。\n");
	return 0;
}
