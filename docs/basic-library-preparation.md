---
title: "Library preparation"
layout: default
parent: Materi dasar
nav_order: 7
---

<p style="text-align: right; font-size: 0.9rem;">
  <a href="https://www.bowo.digital/" style="font-weight: bold;">← Beranda</a>
  &nbsp;&nbsp;|&nbsp;&nbsp;
  <a href="https://www.bowo.digital/docs/part1.html" style="font-weight: bold;">Pilih materi →</a>
</p>

<h1 style="text-align: center; font-size: 2.5rem; font-weight: bold; margin-bottom: 0.5rem;">
  <a href="https://www.bowo.digital/docs/basic-library-preparation.html" style="text-decoration: none; color: inherit;">
    Library preparation
  </a>
</h1>

<p style="text-align: center; font-size: 1.2rem;">
  Oleh <a href="https://www.bowo.digital/docs/bio.html" target="_blank">Agus Wibowo</a>
</p>

<div style="text-align: center; margin-bottom: 1.5rem;">
  <img src="https://cdn.the-scientist.com/assets/articleNo/69901/aImg/45652/lab-scene-epmotion-5073m-the-gripper-tower-an-cmyk-jpg-l.jpg" alt="Library preparation" style="max-width: 100%; height: auto;">
</div>

<br>

# Daftar isi

-   [Pendahuluan](#pendahuluan)
-   [Prinsip dasar library preparation]()
-   [Tahapan utama library preparation]()
-   [Jenis-jenis library preparation]()
-   [Faktor yang mempengaruhi library preparation]()
-   [Library preparation untuk berbagai jenis sekuensing]()
-   [Otomatisasi dalam library preparation]()
-   [Perkembangan terkini]()
-   [Tantangan dan solusi dalam library preparation]()
-   [Panduan praktis memilih metode]()
-   [Kesimpulan]()

# Pendahuluan

Dalam beberapa tahun terakhir, perkembangan teknologi sekuensing terutama *Next-Generation Sequencing* (NGS) telah merevolusi cara kita memahami genom organisme dan menganalisis materi genetik. Di antara tahapan kunci dalam proses ini, *library preparation* (persiapan perpustakaan) memegang peran sentral karena menentukan kualitas hasil sekuensing. Proses ini melibatkan modifikasi DNA/RNA melalui serangkaian langkah biokimia untuk menghasilkan fragmen genetik terorganisir yang siap dibaca oleh mesin sekuensing, mirip seperti buku yang disusun rapi di perpustakaan.

*Library preparation* menjadi penghubung vital antara sampel biologis mentah dan data genetik yang dapat dianalisis. Tingkat keberhasilan proyek sekuensing sangat bergantung pada kualitas library, sehingga pemahaman mendalam tentang tahapan ini penting bagi peneliti genomik. Artikel ini akan mengulas prinsip dasar, tahapan teknis, variasi metode, serta inovasi terkini dalam *library preparation*. Pembahasan juga mencakup faktor penentu kualitas library dan strategi optimasi untuk memastikan hasil sekuensing yang akurat dan efisien.

# Prinsip dasar *library preparation*

Library preparation pada dasarnya adalah proses mengubah asam nukleat (DNA atau RNA) menjadi bentuk yang kompatibel dengan platform sekuensing. Tujuan utamanya adalah menciptakan kumpulan fragmen asam nukleat yang mewakili sampel asli dengan akurat, sekaligus memiliki struktur yang bisa dikenali dan diproses oleh mesin sekuensing.

## Mengapa *library preparation* diperlukan?

Teknologi sekuensing modern tidak dapat langsung membaca molekul DNA/RNA dalam ukuran penuh. Sebagai contoh, genome manusia terdiri dari sekitar 3 miliar pasang basa, sementara kemampuan teknologi sekuensing saat ini hanya mampu membaca beberapa ratus hingga beberapa ribu basa dalam satu kali proses (tergantung platform yang digunakan). Oleh karena itu, materi genetik perlu dipecah menjadi fragmen-fragmen yang lebih kecil.

Selain itu, mesin sekuensing membutuhkan struktur khusus pada ujung-ujung fragmen DNA/RNA, seperti adapter (penghubung) dan barcode (penanda). Adapter berfungsi untuk menempelkan fragmen DNA ke permukaan *flow cell* pada mesin sekuensing, sementara barcode memungkinkan identifikasi sampel ketika beberapa library digabungkan dalam satu run sekuensing (*multiplexing*).

<figure style="text-align: center;">
  <img src="https://nextgen.mgh.harvard.edu/images/LibStructure.png" alt="Illumina Library Structure" style="width: 85%;">
  <figcaption style="font-size: 0.95em; margin-top: 8px; text-align: left;">
    <strong>Struktur library DNA pada platform Illumina.</strong> Setiap fragmen DNA target memiliki adaptor khusus di kedua ujungnya: adaptor P5 dan P7. Adaptor ini memungkinkan fragmen menempel ke flow cell dan dikenali oleh mesin sekuensing. Selain itu, ada elemen-elemen penting lainnya, seperti primer Rd1 dan Rd2 yang digunakan untuk inisiasi proses sekuensing, Index (barcoding) untuk mengidentifikasi sampel berbeda dalam satu run sekuensing (multiplexing), dan primer index seq untuk membaca bagian index selama sekuensing.
    <br><em>Sumber gambar:</em> <a href="https://nextgen.mgh.harvard.edu/" target="_blank">Nextgen.mgh.harvard.edu</a>
  </figcaption>
</figure>

Kualitas *library* sangat menentukan kualitas data sekuensing yang dihasilkan. *Library preparation* yang baik akan menghasilkan fragmen dengan ukuran yang sesuai, distribusi yang merata, dan *coverage* yang cukup di seluruh genom. Sebaliknya, *library preparation* yang buruk dapat menyebabkan bias dalam representasi beberapa region genom, coverage yang tidak merata, atau bahkan hilangnya informasi dari region genom tertentu.

# Tahapan utama library preparation

Meskipun detail proses library preparation bervariasi tergantung pada jenis sampel dan platform sekuensing yang digunakan, secara umum terdapat 7 tahapannya yaitu:

## Ekstraksi dan purifikasi asam nukleat

Tahap pertama adalah mendapatkan DNA atau RNA berkualitas tinggi dari sampel biologis. Ini melibatkan pemecahan sel (lisis), pemisahan asam nukleat dari komponen sel lainnya, dan pemurnian untuk menghilangkan kontaminan. Kualitas materi genetik awal sangat memengaruhi keberhasilan tahapan selanjutnya.

## Fragmentasi DNA/RNA

Pada tahap ini, DNA/RNA dipecah menjadi potongan-potongan yang lebih kecil dengan ukuran yang sesuai untuk platform sekuensing yang digunakan. Fragmentasi dapat dilakukan dengan beberapa metode:

-   **Fragmentasi fisik**: Menggunakan sonikasi, nebulisasi, atau *shearing* hidrodinamik.
-   **Fragmentasi enzimatik**: Menggunakan enzim seperti *DNase I* atau *endonuklease* pembatas lainnya.
-   **Fragmentasi kimia**: Menggunakan panas dan kation divalen.

Misalnya, pada protokol persiapan library Illumina, fragmentasi biasanya menghasilkan potongan DNA berukuran 200-800bp yang merupakan ukuran ideal untuk sekuensing jangka pendek (*short-read sequencing*).

## End repair dan dA-tailing

Fragmentasi sering kali menghasilkan ujung DNA yang tidak rata, sehingga perbaikan (*end repairing*) pada ujung-ujung DNA ini diperlukan. Tahap *end repair* melibatkan proses pemotongan dan pengisian untuk menghasilkan ujung tumpul (blunt ends). Selanjutnya, proses dA-tailing dikerjakan untuk menambahkan nukleotida adenin (A) ke ujung 3' fragmen, yang nantinya akan berpasangan dengan nukleotida timin (T) pada adapter.

## Ligasi adapter

Adapter, yang merupakan oligonukleotida pendek dengan sekuens khusus, ditambahkan ke kedua ujung fragmen DNA. Adapter ini berfungsi sebagai:

1.  Tempat menempelnya primer untuk amplifikasi PCR
2.  Tempat pengikatan fragmen ke flow cell
3.  Tempat inisiasi reaksi sekuensing

## Seleksi ukuran

Tahap ini memastikan bahwa fragmen yang akan disekuensing memiliki ukuran yang seragam dan optimal. Seleksi ukuran biasanya dilakukan menggunakan elektroforesis gel, AMPure XP beads, atau metode lainnya. Misalnya, protokol Illumina DNA Prep biasanya menyertakan tahap purifikasi menggunakan AMPure beads untuk memilih fragmen dalam rentang ukuran tertentu.

## Amplifikasi PCR (opsional)

Untuk meningkatkan jumlah materi yang tersedia untuk sekuensing, *library* sering diamplifikasi menggunakan PCR. Tahap ini juga dapat digunakan untuk menambahkan barcode atau indeks ke *library*, yang memungkinkan *multiplexing* sampel dalam satu run sekuensing.

## Validasi *Library*

Sebelum sekuensing, kualitas dan kuantitas *library* perlu diperiksa menggunakan metode seperti elektroforesis kapiler (Bioanalyzer/TapeStation) dan kuantifikasi fluorometrik (Qubit). Validasi ini memastikan bahwa *library* memiliki konsentrasi yang cukup, ukuran fragmen yang sesuai, dan tidak mengandung kontaminan seperti dimer adapter.

Rangkuman proses-proses penting dari *library preparation* diilustrasikan pada gambar berikut.

<figure style="text-align: center;">
  <img src="images/lib-prep.png" alt="Illumina Library prep" style="width: 85%;">
  <figcaption style="font-size: 0.95em; margin-top: 8px; text-align: left;">
    <strong>Ilustrasi proses library preparation.</strong>.
    <br><em>Sumber gambar:</em> <a href="https://link.springer.com/book/10.1007/978-3-030-62490-3" target="_blank">NGS and Data Analysis</a>
  </figcaption>
</figure>

# 

--- Sekian ---

# Referensi

-   [Melanie Kappelmann-Fenzl. (2021). NGS and Data Analysis. Springer Cham](https://doi.org/10.1007/978-3-030-62490-3)