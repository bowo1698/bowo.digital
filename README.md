---
layout: default
title: Home
nav_order: 1
description: "Panduan belajar bioinformatika dengan fokus kajian pada akuakultur."
permalink: /
---

<h1 style="text-align: center; font-size: 2.5rem; font-weight: bold; margin-bottom: 0.5rem;">
  <a href="https://www.bowo.digital/" style="text-decoration: none; color: inherit;">
    Bioinformatic Notes
  </a>
</h1>

<p style="text-align: center; font-size: 1.2rem;">
  Oleh <a href="https://www.bowo.digital/docs/bio.html" target="_blank">Agus Wibowo</a>
</p>

<div style="text-align: center; margin-bottom: 1.5rem;">
  <img src="https://www.genomicsengland.co.uk/assets/imager/images/Technology/315155/Bioinformatics-and-data-hands_6c0c164bd2b597ee32b68b8b5755bd2e.jpg" alt="Tools instalation" style="max-width: 100%; height: auto;">
</div>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![GitHub issues](https://img.shields.io/github/issues/bowo1698/bowo.digital) ![GitHub repo size](https://img.shields.io/github/repo-size/bowo1698/bowo.digital) [![DOI](https://zenodo.org/badge/243280413.svg)](https://zenodo.org/badge/latestdoi/243280413)

## Selamat datang 👋

### 🔬 Sekilas tentang bioinformatika

Bioinformatika adalah ilmu interdisipliner yang mengintegrasikan prinsip biologi, statistika, matematika, dan komputasi untuk menjawab tantangan biologi modern. Ketertarikan saya bermula dari perkembangan ilmu akuakultur yang kini tak lagi terbatas pada pemeliharaan ikan di kolam-kolam ataupun analisis di laboratorium basah (*wet lab*), tetapi merambah ke analisis data molekuler berskala besar. Di sinilah peran krusial bioinformatika, mengolah, menganalisis, dan memodelkan data biologis untuk mengungkap mekanisme kehidupan di tingkat molekuler. Ini seperti mempelajari suatu materi fisik, tetapi di level kuantum!

<div style="text-align: center;">
  <img src="assets/fish-with-com.jpg" alt="Fish With Code" style="width: 90%; height: auto;">
</div>

Yang menarik adalah, karena kita bermain dengan "ikan", eksplorasi biologi menjadi lebih *flexible* dan membuka banyak ruang inovasi. Misalnya, kita bisa menelusuri gen-gen yang terkait dengan pertumbuhan cepat, ketahanan terhadap penyakit, atau efisiensi pakan, semuanya hanya melalui data DNA. Bioinformatika memberi kita kacamata baru untuk melihat ikan bukan sekadar objek konsumsi, tetapi sebagai kumpulan informasi genetik yang bisa kita manfaatkan untuk perbaikan populasi secara presisi dan berkelanjutan, hingga akhirnya bisa menjadi model atau landasan penerapan untuk level yang lebih tinggi, seperti manusia.

Perkembangan mesin sekuensing generasi mutakhir seperti [Illumina](https://en.wikipedia.org/wiki/Illumina,_Inc.) (*short-read*) dan [PacBio](https://en.wikipedia.org/wiki/Pacific_Biosciences) (*long-read*) telah menghasilkan data biologis dalam volume yang luar biasa besar. Bioinformatika menjadi kunci untuk mengekstrak wawasan dari kompleksitas ini, terutama dalam [tiga bidang utama terkait dengan akuakultur](https://doi.org/10.1007/978-981-97-8553-7_11):

-   ***Genomics***: Bioinformatika digunakan untuk mempercepat program pemuliaan berbasis genom melalui identifikasi variasi berbagai marker polimorfik seperti [mikrosatelit](https://en.wikipedia.org/wiki/Microsatellite), [SNP](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism), atau bahkan yang terbaru, [mikrohaplotipe](https://pubmed.ncbi.nlm.nih.gov/30347322/). Yang paling dikejar adalah merancang *guide RNA* pada teknologi **CRISPR/Cas9** guna mengedit gen target (jika suatu sifat dikendalikan oleh sedikit gen, *non-polygenic*), seperti meningkatkan toleransi oksigen rendah pada udang atau ketahanan terhadap virus di ikan kerapu.

-   ***Metagenomics***: Bioinformatika dapat membantu menganalisis keragaman mikroba dalam ekosistem akuatik (misalnya, air tambak atau saluran pencernaan ikan) untuk mengidentifikasi komunitas bakteri menguntungkan. Dari data ini, kita dapat mendesain probiotik spesifik yang dapat meningkatkan kesehatan ikan, mengurangi ketergantungan pada antibiotik, dan mengoptimalkan kualitas air.

-   ***Transcriptomics***: Bioinformatika dapat memetakan ekspresi gen organisme akuatik di bawah kondisi tertentu (stres, penyakit, infeksi, perlakuan khusus, dll) untuk merancang program pemuliaan berbasis marker genetik. Contohnya, mengembangkan strain ikan dengan pertumbuhan cepat atau ketahanan terhadap patogen seperti *Vibrio*.

### 🎓 Tujuan website

Website ini dirancang sebagai wadah kolaborasi, di mana siapa pun dapat berkontribusi dalam berbagi pengetahuan seputar bioinformatika dan akuakultur. Tujuannya sederhana, membangun sumber belajar yang terbuka, mudah diakses, dan relevan dengan kebutuhan praktis, khususnya dalam konteks akuakultur di Indonesia yang memiliki potensi luar biasa.

Seluruh konten, termasuk tutorial dan panduan, disajikan dalam bahasa Indonesia dengan gaya penyampaian yang sederhana agar mudah dipahami oleh berbagai kalangan. Semua materi dapat dipelajari secara gratis, menggunakan data *open access* dan perangkat lunak *open source*, sehingga tidak ada hambatan biaya maupun teknis untuk ikut serta.

Dengan semangat kolaboratif, platform ini mengajak semua pihak untuk saling berbagi dan belajar, demi mewujudkan akuakultur yang berkelanjutan dan berbasis ilmu pengetahuan di Indonesia.

> "*We are a great nation - what we need are people willing to make a difference, even from behind a screen full of code.*"

## 💻 Catatan teknis

Karena sebagian besar tools bioinformatika dikembangkan dan dioptimalkan untuk berjalan di lingkungan UNIX terutama [LINUX](https://en.wikipedia.org/wiki/Linux) dan MACOS, saya sangat menyarankan rekan-rekan pengguna Windows untuk mempertimbangkan "bermigrasi" ke sistem operasi berbasis Linux ini, seperti [Ubuntu](https://ubuntu.com/), [Debian](https://www.debian.org/), atau [Fedora](https://fedoraproject.org/).

Mengapa demikian?

1.  **Kompatibilitas** – Tools seperti HISAT2, STAR, Salmon, dan RSEM berjalan lebih stabil dan langsung didukung di lingkungan UNIX. "Sangat mudah" dalam hal instalasi dan eksekusinya

2.  **Performa** – Linux cenderung menggunakan sumber daya sistem lebih efisien dibandingkan Windows, terutama dengan adanya sistem *swap* membuatnya ideal untuk komputasi bioinformatika, bahkan di perangkat berspesifikasi rendah.

3.  **Minim gangguan** – Sistem Linux jauh lebih tahan (bahkan hampir mustahil) terhadap malware, dan mendukung scripting otomatis yang efisien (Bash/Shell).

4.  **Ekosistem terbuka** – Hampir semua perangkat lunak bioinformatika bersifat *open-source* (gratis) dan dikembangkan di lingkungan Linux, membuat dokumentasi dan komunitas dukungannya jauh lebih luas.

<div style="text-align: center;">
  <img src="assets/linux-win-mac.jpeg" alt="Fish With Code" style="width: 90%; height: auto;">
</div>

> Seluruh materi dalam website ini disusun dengan asumsi bahwa Anda menjalankan Linux secara langsung (*full system install*), bukan melalui Windows. Panduan penggunaan *Windows Subsystem for Linux* (WSL) hanya disediakan sebagai opsi tambahan untuk keperluan eksplorasi atau pembelajaran awal. Namun, dalam praktiknya, WSL memiliki banyak keterbatasan dalam hal kompatibilitas dan kinerja, sehingga tidak direkomendasikan untuk analisis bioinformatika yang sesungguhnya.

Untuk perangkat keras, saya merancang semua tutorial agar dapat dijalankan pada komputer/laptop dengan spesifikasi medium (*laptop dana pelajar* 😁), selama menggunakan sistem operasi berbasis Linux [cek disini](docs/basic-linux-in-win.md). Hal ini memungkinkan siapapun tetap dapat mengikuti alur tutorial ini, mulai dari unduhan data, pra-pemrosesan, analisis, hingga visualisasi, tanpa memerlukan perangkat mahal!!! (meskipun lebih baik jika setiap kampus punya *High Performance Computing* - HPC, yuk bikin!)

## 🚀 Materi

Untuk mulai belajar, silahkan pilih tutorial mana yang anda inginkan.

### [1. Materi dasar](docs/part1.md)

-   [Linux di dalam Windows](docs/basic-linux-in-win.md)
-   [Memahami hierarki Linux](docs/basic-hierarki-linux.md)
-   [Instalasi tools bioinformatika](docs/basic-instalasi-tools.md)
-   [Konsep-konsep dasar biologi molekuler](docs/basic-kosep-biomol.md)
-   [Teknik dasar biologi molekuler](docs/basic-teknik-biomol.md)
-   [Tren dan perkembangan aplikasi biologi molekuler dan bioinformatika untuk akuakultur](docs/basic-tren-terbaru.md)
<!--   [Marker-marker genetik](docs/basic-marker-genetik.md) --->
<!--   [Pengenalan dasar-dasar Bash](docs/basic-bash.md) --->
<!--   [Pengenalan dasar-dasar R](docs/basic-R.md) --->
<!--   [Pengenalan dasar-dasar python](docs/basic-python.md) --->
<!--   [Dokumentasi melalui Markdown](docs/basic-markdown.md) --->
<!--   [Eksplorasi database genomic](docs/basic-eksplorasi-database.md) --->
<!--   [Cara download file *sequencing*](docs/basic-curl-wget.md) --->
<!--   [Format file yang digunakan dalam bioinformatika](docs/basic-file_format.md) --->
<!--   [Cara kerja *sequencing*](docs/basic-sequencing.md) --->
<!--   [Mengenal tentang studi omics](docs/basic-omics.md) --->
<!--   [Kosep dasar *alignment*](docs/basic-alignment.md) --->
<!--   [Statistik untuk biologi modern](docs/basic-stats.md) --->

### [2. Materi lanjutan 1](docs/part2.md)

<!--   [Quality control hasil *sequencing*](docs/tutorial-quality-control.md) --->
<!--   [*Alignment*](docs/tutorial-alignment.md) --->
<!--   [Pipeline metagenomics](docs/tutorial-metagenomics.md) --->
<!--   [Pipeline RNA-Seq](docs/tutorial-RNA-Seq.md) --->
<!--   [Pipeline genomic varians - SNP discovery](docs/tutorial-SNP-discovery.md) --->
<!--   [Pipeline genomic varians - haplotype block discovery](docs/tutorial-haplotype.md) --->
<!--   [Pipeline genomic varians - microhaplotype discovery](docs/tutorial-microhaplotype.md) --->

### [3. Materi lanjutan 2](docs/part3.md)

<!--   [Genetika populasi](docs/adv-GWAS.md) --->
<!--   [Pemuliaan berbasis genomic - GWAS](docs/adv-GWAS.md) --->
<!--   [Pemuliaan berbasis genomic - estimasi nilai breeding](docs/adv-EBV.md) --->

<!-- komentar -->
<p style="text-align: right; font-size: 0.9rem;">
  <a href="https://www.bowo.digital/" style="font-weight: bold;">← Beranda</a>
  &nbsp;&nbsp;|&nbsp;&nbsp;
  <a href="https://www.bowo.digital/docs/basic-tren-terbaru.html" style="font-weight: bold;">Konten berikutnya →</a>
</p>

