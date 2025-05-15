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


<div style="text-align: center;">
  <h2 style="font-size: 2rem;">Selamat datang 👋</h2>
</div>

<div style="text-align: justify;">
  <h3>🔬 Sekilas tentang bioinformatika</h3>
  <p>
    Bioinformatika adalah ilmu interdisipliner yang mengintegrasikan prinsip biologi, statistika, matematika, dan komputasi untuk menjawab tantangan biologi modern. Ketertarikan saya bermula dari perkembangan ilmu akuakultur yang kini tak lagi terbatas pada pemeliharaan ikan di kolam-kolam ataupun analisis di laboratorium basah (<em>wet lab</em>), tetapi merambah ke analisis data molekuler berskala besar. Di sinilah peran krusial bioinformatika, mengolah, menganalisis, dan memodelkan data biologis untuk mengungkap mekanisme kehidupan di tingkat molekuler. Ini seperti mempelajari suatu materi fisik, tetapi di level kuantum!
  </p>

  <div style="text-align: center;">
    <img src="assets/fish-with-com.jpg" alt="Fish With Code" style="width: 90%; height: auto;">
  </div>

  <p>
    Yang menarik adalah, karena kita bermain dengan "ikan", eksplorasi biologi menjadi lebih <em>flexible</em> dan membuka banyak ruang inovasi. Misalnya, kita bisa menelusuri gen-gen yang terkait dengan pertumbuhan cepat, ketahanan terhadap penyakit, atau efisiensi pakan, semuanya hanya melalui data DNA. Bioinformatika memberi kita kacamata baru untuk melihat ikan bukan sekadar objek konsumsi, tetapi sebagai kumpulan informasi genetik yang bisa kita manfaatkan untuk perbaikan populasi secara presisi dan berkelanjutan, hingga akhirnya bisa menjadi model atau landasan penerapan untuk level yang lebih tinggi, seperti manusia.
  </p>

  <p>
    Perkembangan mesin sekuensing generasi mutakhir seperti <a href="https://en.wikipedia.org/wiki/Illumina,_Inc.">Illumina</a> (<em>short-read</em>) dan <a href="https://en.wikipedia.org/wiki/Pacific_Biosciences">PacBio</a> (<em>long-read</em>) telah menghasilkan data biologis dalam volume yang luar biasa besar. Bioinformatika menjadi kunci untuk mengekstrak wawasan dari kompleksitas ini, terutama dalam <a href="https://doi.org/10.1007/978-981-97-8553-7_11">tiga bidang utama terkait dengan akuakultur</a>:
  </p>

  <ul>
    <li>
      <strong><em>Genomics</em></strong>: digunakan untuk mempercepat program pemuliaan berbasis genom melalui identifikasi marker seperti <a href="https://en.wikipedia.org/wiki/Microsatellite">mikrosatelit</a>, <a href="https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism">SNP</a>, atau <a href="https://pubmed.ncbi.nlm.nih.gov/30347322/">mikrohaplotipe</a>. Bioinformatika juga membantu merancang <em>guide RNA</em> untuk teknologi <strong>CRISPR/Cas9</strong> dalam mengedit gen target.
    </li>
    <li>
      <strong><em>Metagenomics</em></strong>: untuk menganalisis komunitas mikroba dalam ekosistem akuatik dan merancang probiotik yang meningkatkan kesehatan ikan dan kualitas air.
    </li>
    <li>
      <strong><em>Transcriptomics</em></strong>: memetakan ekspresi gen dalam kondisi tertentu untuk mendukung pemuliaan berbasis marker genetik, seperti ketahanan terhadap patogen.
    </li>
  </ul>

  <h3>🎓 Tujuan website</h3>
  <p>
    Website ini dirancang sebagai wadah kolaborasi untuk berbagi pengetahuan seputar bioinformatika dan akuakultur. Tujuannya adalah membangun sumber belajar yang terbuka, mudah diakses, dan relevan secara praktis, khususnya untuk konteks akuakultur Indonesia.
  </p>
  <p>
    Seluruh konten disajikan dalam bahasa Indonesia dengan penyampaian yang sederhana agar dapat dipahami oleh berbagai kalangan. Semua materi gratis, menggunakan data <em>open access</em> dan perangkat lunak <em>open source</em>, tanpa hambatan biaya atau teknis.
  </p>
  <p>
    Dengan semangat kolaboratif, platform ini mengajak semua pihak untuk saling berbagi dan belajar demi mewujudkan akuakultur yang berkelanjutan dan berbasis sains di Indonesia.
  </p>

  <blockquote>
    <em>"We are a great nation - what we need are people willing to make a difference, even from behind a screen full of code."</em>
  </blockquote>

  <h2>💻 Catatan teknis</h2>
  <p>
    Karena sebagian besar tools bioinformatika dikembangkan untuk lingkungan UNIX, terutama <a href="https://en.wikipedia.org/wiki/Linux">Linux</a> dan macOS, saya menyarankan pengguna Windows mempertimbangkan sistem operasi Linux seperti <a href="https://ubuntu.com/">Ubuntu</a>, <a href="https://www.debian.org/">Debian</a>, atau <a href="https://fedoraproject.org/">Fedora</a>.
  </p>
  <ol>
    <li><strong>Kompatibilitas</strong> – tools seperti HISAT2, STAR, Salmon, dan RSEM lebih stabil di UNIX.</li>
    <li><strong>Performa</strong> – Linux efisien dalam penggunaan sumber daya sistem.</li>
    <li><strong>Minim gangguan</strong> – aman dari malware dan mendukung scripting otomatis.</li>
    <li><strong>Ekosistem terbuka</strong> – dokumentasi dan dukungan komunitas luas.</li>
  </ol>

  <div style="text-align: center;">
    <img src="assets/linux-win-mac.jpeg" alt="Linux vs Windows vs Mac" style="width: 90%; height: auto;">
  </div>

  <blockquote>
    Seluruh materi di website ini diasumsikan dijalankan di sistem Linux penuh, bukan WSL. Panduan WSL hanya untuk eksplorasi, dan tidak direkomendasikan untuk analisis nyata karena keterbatasan kompatibilitas dan performa.
  </blockquote>

  <p>
    Tutorial ini dirancang agar dapat dijalankan di laptop spesifikasi menengah selama menggunakan Linux. Lihat <a href="docs/basic-linux-in-win.md">di sini</a> untuk panduan sistem.
  </p>

  <h2>🚀 Materi</h2>
  <p>
    Untuk mulai belajar, silakan pilih tutorial mana yang Anda inginkan:
  </p>
  <h3><a href="docs/part1.md">1. Materi dasar</a></h3>

  <hr style="margin: 2rem 0;">
  <p style="text-align: center; font-size: 0.95rem;">
    <a href="/" style="font-weight: bold;">← Beranda</a>
    &nbsp;&nbsp;|&nbsp;&nbsp;
    <a href="docs/part2.md" style="font-weight: bold;">Lanjut ke Materi Berikutnya →</a>
  </p>
</div>

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


