---
title: "Quality control hasil sequencing"
layout: default
parent: Skill Teknis
nav_order: 2
---

<p style="text-align: right; font-size: 0.9rem;">
  <a href="https://www.bowo.digital/" style="font-weight: bold;">← Beranda</a>
  &nbsp;&nbsp;|&nbsp;&nbsp;
  <a href="https://www.bowo.digital/docs/part1.html" style="font-weight: bold;">Pilih materi →</a>
</p>

<h1 style="text-align: center; font-size: 2.5rem; font-weight: bold; margin-bottom: 0.5rem;">
  <a href="https://www.bowo.digital/docs/teknis-quality-control.html" style="text-decoration: none; color: inherit;">
    Quality control hasil sequencing
  </a>
</h1>

<p style="text-align: center; font-size: 1.2rem;">
  Oleh <a href="https://www.bowo.digital/docs/bio.html" target="_blank">Agus Wibowo</a>
</p>

<div style="text-align: center; margin-bottom: 1.5rem;">
  <img src="https://www.integra-biosciences.com/sites/default/files/styles/large/public/images/data-analysis.jpg?itok=DX2veo4F" alt="QC" style="max-width: 100%; height: auto;">
</div>

<br>

# Daftar isi

-   [Apa itu quality control dan mengapa penting?](#apa-itu-quality-control-dan-mengapa-penting)
-   [Mendownload data sekuen]()
-   [FASTQC]()
-   [SAM & BAM: Dari read ke posisi di genom]()
-   [GFF/GTF & BED: Menyematkan makna pada genom]()
-   [VCF: Menyimpan informasi variasi genetik]()
-   [SRA: Gudangnya data NGS dunia]()
-   [Memahami format = Menguasai analisis]()

>   **Catatan**: Dalam tutorial ini, kita akan bekerja menggunakan pemrograman Bash dalam lingkungan Linux. Anda dapat menggunakan VS code baik untuk membuat script R maupun menjalankan script Bash tersebut via terminal.

# Apa itu *quality control* dan mengapa penting?

Bayangkan Anda mengambil foto dengan kamera harga ratusan juta, tapi lensa yang digunakan berembun. Hasilnya? Gambar blur yang nyaris tak bisa dilihat. Begitu juga dengan data sekuensing DNA yang kualitasnya buruk—tak peduli seberapa canggih analisisnya, hasilnya tetap tidak akan optimal.

*Next-Generation Sequencing* (NGS) telah membuka jendela baru untuk mengintip rahasia genetik makhluk hidup. Dalam sekali eksperimen, jutaan fragmen DNA berhasil dibaca! Tapi jumlah besar ini tidak otomatis menjamin kualitas. Seperti memilah berlian dari kerikil, *quality contro*l (QC) menjadi langkah krusial untuk memastikan data yang kita analisis benar-benar berharga.

# Mendownload data sekuen

Dalam tutorial ini, kita menggunakan data dari studi oleh [Robledo et al. (2014)](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-1149) yang meneliti respons transkriptomik ikan turbot (*Scophthalmus maximus*) terhadap infeksi *Enteromyxum scophthalmi*, penyebab enteromikosis. Penelitian ini menganalisis perubahan ekspresi gen di ginjal, limpa, dan usus menggunakan RNA-seq berbasis Illumina HiSeq 2000. 

Dataset lengkap dapat diakses melalui [ENA (European Nucleotide Archive)](https://www.ebi.ac.uk/ena/browser/view/PRJNA269386) dengan kode proyek PRJNA269386. Karena ukuran total data mencapai sekitar 27 GB, dalam tutorial ini kita hanya menggunakan satu file sequencing yaitu SRR1695153, yang berasal dari organ limpa. File ini telah dimodifikasi untuk keperluan edukasi sehingga ukurannya menjadi lebih kecil (56 MB) dan lebih mudah diolah. Bagi yang tertarik mempelajari proses manipulasi data tersebut, kode lengkapnya dapat didownload <a href="/code/fastq_manip.py" download>di sini</a>.

Buat forlder baru untuk menyimpan dan mendownload file FASTQ yang sudah dikompresi dengan gunzip dengan cara:

```bash
# bash

# buat folder baru
mkdir -p data_fastq/

# masuk ke direktori data_fastq/
cd data_fastq/

# download data SRR
wget https://huggingface.co/datasets/Antijokowisme16/fastqc_data/resolve/main/SRR1695153_R1_edit.fastq.gz
wget https://huggingface.co/datasets/Antijokowisme16/fastqc_data/resolve/main/SRR1695153_R2_edit.fastq.gz

# ubah akses ke read-only
chmod u-w *fastq.gz
```

> **Catatan:** Setiap hasil sequencing yang menggunakan *paired-end reads* maka akan selalu memiliki 2 file reads (R1 dan R2) untuk setiap sampel. 

# FASTQC

Setalah data didownload dan pada terminal kita berada di direktori `data_fastq/`, maka langkah selanjutnya adalah menginstall tools [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). FASTQC adalah tool yang digunakan untuk melakukan *quality control* awal terhadap data *sequencing*, dengan memberikan ringkasan visual mengenai kualitas base, *adapter contamination*, *overrepresented sequences*, dan metrik penting lainnya dalam file FASTQ.

Install FASTQC melalui anaconda Python dengan langkah-langkah sebagai berikut:

```bash
# bash
# buat env khusus untuk QC data sequencing dengan nama misalnya: qc_sequencing
conda create --name qc_sequence python=3.11

# aktifkan env qc_sequencing
conda activate qc_sequencing

# install fastqc
conda install fastqc
```

Untuk melakukan quality control dengan FASTQC pada data FASTQ, cukup gunakan perintah berikut setelah membuat direktori baru untuk menyimpan file hasil QC.

```bash
# bash
# buat direktori hasil QC
mkdir -p fastqc

# jalankan fastqc
fastqc *.fastq.gz -o fastqc/
```

> **Catatan:** Setiap perintah yang menggunakan simbol `*` diikuti dengan jenis file, misalnya `*.fastq.gz` maka perintah tersebut akan memproses semua file dengan format .fastq.gz. Teknik ini sangat bermanfaat jika kita ingin memproses banyak file dengan jenis file yang sama.

Setelah menjalankan FASTQC, di dalam folder `fastqc` akan terdapat dua jenis file: file berakhiran `.html` yang berisi laporan visual hasil analisis kualitas, dan file `.zip` yang berisi data mentah serta detail dari laporan tersebut. Untuk keperluan evaluasi, kita cukup membuka file `.html` menggunakan browser seperti Chrome atau Safari, maka akan tambil hasil seperti ini.

<div style="text-align: center;">
  <img src="images/fastqc.png" alt="Contoh hasil FASTQC" style="width: 90%; height: auto;">
</div>

Anda juga bisa observasi langsung hasil FASTQC tanpa harus menjalankan perintah `fastqc` dengan klik link berikut: [SRR1695153_R1](results/SRR1695153_1_edit_fastqc.html) dan [SRR1695153_R2](results/SRR1695153_2_edit_fastqc.html)

## Interpretasi

