# fusion-metrics (Tiếng Việt)

[English](README.md) · [中文](README.zh.md)

**Bộ chỉ số đánh giá chất lượng dung hợp ảnh (image fusion), dịch từng dòng từ MATLAB và kiểm chứng khớp đến sai số làm tròn.**

Tác giả: **[Vũ Việt Anh](https://vuvietanh2309.github.io/)**, AI Engineer, học viên Thạc sĩ Khoa học Máy tính, Đại học Bách khoa Hà Nội (HUST)

## Vì sao cần thư viện này

Các bài báo về image fusion thường báo cáo Q<sup>AB/F</sup>, MI, NCIE, FMI, VIF… nhưng mỗi nhóm tính bằng một bộ script
MATLAB hoặc một bản Python tự viết, và các bản này lệch nhau mà không ai để ý. Nguyên nhân thường là hệ số đổi ảnh xám
khác nhau, `conv2(...,'same')` lệch một pixel, `uint8` bị bão hoà trong histogram, hay quy tắc làm mảnh cạnh khác nhau.
Kết quả là cùng một phương pháp nhưng mỗi bài ra một số.

`fusion-metrics` cung cấp **một bản cài đặt chuẩn** mà ai cũng `pip install` được:

- **Trung thành với MATLAB.** Mỗi hàm dịch từng dòng từ một file `.m` tham chiếu, giữ nguyên tên hàm, thứ tự tham số và
  các chi tiết số học.
- **Đã kiểm chứng.** Self-test so từng giá trị với kết quả của chính file `.m`; độ lệch ≤ 1e-12. Test chạy tự động trên
  Python 3.10–3.13.
- **Dùng được ngay cho benchmark.** Gọi từ Python cho một ảnh, hoặc chấm cả thư mục kết quả của một phương pháp bằng dòng
  lệnh và xuất CSV.

## Cài đặt và dùng thử

```bash
pip install fusion-metrics
```

```python
import fusion_metrics as fm

A, B, F = fm.read_gray("MRI.png"), fm.read_gray("PET.png"), fm.read_gray("fused.png")
fm.Q_ABF(A, B, F)                     # một metric
fm.evaluate(A, B, F, "information")   # một nhóm metric
fm.evaluate(A, B, F)                  # cả 35 metric
```

Dòng lệnh:

```bash
fusion-metrics MRI.png PET.png fused.png --preset feature
fusion-metrics --dirs data/MRI data/PET results/MyMethod --csv mymethod.csv   # chấm cả thư mục
python verify/selftest.py                                                     # tự kiểm tra, in "ALL OK"
```

Các nhóm metric có sẵn: `image`, `information`, `feature`, `structural`, `perception`, `similarity`, `ivif`, `liu2012` (đủ 12 metric của Liu et al., TPAMI 2012), `all`. Danh sách đầy đủ, nguồn trích dẫn của
từng metric và các quy ước làm thay đổi kết quả (đổi ảnh xám, kiểu dữ liệu `uint8`, `edge()`…) có trong
[README tiếng Anh](README.md).

## Lưu ý quan trọng

- Đọc ảnh bằng `read_gray` (giống MATLAB `imread` + `rgb2gray`) để ra đúng số của MATLAB. `cv2.imread(path, 0)` cho kết
  quả khác.
- `Q_AG` khác `AG_VIFB`, và `Q_VIFF` (Han et al.) khác `VIF_Sheikh` (Sheikh & Bovik). Khi báo cáo, ghi rõ dùng metric nào.
- `MS_SSIM` thực ra là SSIM một thang (lấy trung bình SSIM(A,F) và SSIM(B,F)), đúng như file `.m` gốc.
- `Q_UIQI` giữ đúng hành vi của `UIQI.m`: file này không đổi ảnh sang double, nên với ảnh `uint8` các phép bình phương bị bão hoà ở 255 như trong MATLAB. Muốn có UIQI theo định nghĩa gốc thì truyền ảnh `float64`.
- `Q_C`, `Q_Y` (SSIM với C = 2e-16) và `Q_CB` (tỉ số của hai ảnh đã lọc FFT) vốn không ổn định số học ngay trong file `.m`: ở các vùng phẳng khác 0 (Q_C, Q_Y) hoặc nơi mẫu số qua 0 (Q_CB), MATLAB, Octave và NumPy cho kết quả khác nhau. Bản Python tính các cửa sổ phẳng bằng số học chính xác cho Q_C/Q_Y; Q_CB chỉ tái lập được đến khoảng 1e-4 giữa các thư viện FFT. Chi tiết ở README tiếng Anh, mục *Verification*.
- Khi dùng, hãy trích dẫn cả thư viện và **bài báo gốc của từng metric**.

## Tác giả

**Vũ Việt Anh**: tác giả toàn bộ thư viện Python, dòng lệnh và hệ thống kiểm chứng.<br/>
Email: [vietanhresearcher@gmail.com](mailto:vietanhresearcher@gmail.com) · ORCID: [0009-0006-6977-5922](https://orcid.org/0009-0006-6977-5922) · Website: [vuvietanh2309.github.io](https://vuvietanh2309.github.io/) · GitHub: [@VuVietAnh2309](https://github.com/VuVietAnh2309)

## Giấy phép

Mã nguồn theo giấy phép MIT. Riêng `Q_FSIM`, `Q_VIFF` và `VIF_Sheikh` là bản dịch từ mã của tác giả khác, chỉ được dùng
cho mục đích học tập và nghiên cứu; xem [NOTICE](NOTICE).

Các file `.m` trong `matlab/` là bản tham chiếu dùng để kiểm chứng; file nào lấy từ toolbox đã công bố thì giữ nguyên phần
ghi chú gốc.

Góp ý, báo lỗi hoặc đóng góp bản port cho các metric còn lại: [GitHub Issues](https://github.com/VuVietAnh2309/fusion-metrics/issues).
