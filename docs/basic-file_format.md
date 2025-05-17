---
title: "Mengenal berbagai format file NGS"
layout: default
parent: Materi dasar
nav_order: 6
---

<p style="text-align: right; font-size: 0.9rem;">
  <a href="https://www.bowo.digital/" style="font-weight: bold;">← Beranda</a>
  &nbsp;&nbsp;|&nbsp;&nbsp;
  <a href="https://www.bowo.digital/docs/part1.html" style="font-weight: bold;">Pilih materi →</a>
</p>

<h1 style="text-align: center; font-size: 2.5rem; font-weight: bold; margin-bottom: 0.5rem;">
  <a href="https://www.bowo.digital/docs/basic-linux-in-win.html" style="text-decoration: none; color: inherit;">
    Mengenal berbagai format file Next Generation Sequencing
  </a>
</h1>

<p style="text-align: center; font-size: 1.2rem;">
  Oleh <a href="https://www.bowo.digital/docs/bio.html" target="_blank">Agus Wibowo</a>
</p>

<div style="text-align: center; margin-bottom: 1.5rem;">
  <img src="https://wallpapers.com/images/hd/linux-pictures-g8nflkx3vjw6v0xq.jpg" alt="Linux in Windows" style="max-width: 100%; height: auto;">
</div>

# Daftar isi

-   [Pendahuluan]()
-   [Apa itu read, mate, dan coverage?]()
-   [FASTA & FASTQ: Format dasar NGS dan kualitasnya]()
-   [SAM & BAM: Dari read ke posisi di genom]()
-   [GFF/GTF & BED: Menyematkan makna pada genom]()
-   [VCF: Menyimpan informasi variasi genetik]()
-   [SRA: Gudangnya data NGS dunia]()
-   [BedGraph & BigWig: Visualisasi dari data kuantitatif]()
-   [Memahami format = Menguasai analisis]()

# Pendahuluan

Bayangkan kamu baru saja mendapatkan data hasil sekuensing dari laboratorium. Kamu pikir: “Wah, sebentar lagi saya bisa tahu gen apa yang rusak atau SNP apa yang unik.” Tapi, begitu buka foldernya… kamu dihadapkan dengan berbagai file berakhiran .fastq, .bam, .gff, .vcf, dan entah apa lagi. “mumet wes” 😵‍💫

Yup, dunia Next-Generation Sequencing (NGS) memang penuh dengan format file yang bagi pemula bisa sangat membingungkan. Tapi jangan khawatir—setiap format punya fungsinya sendiri. Seperti alat-alat di toolbox, mereka bekerja bersama untuk menjawab pertanyaan besar di biologi molekuler: variasi genetik, ekspresi gen, struktur genom, dan banyak lagi.

Jadi, konten ini akan membawa kamu berkenalan dengan “bahasa” NGS: format file yang menjadi jembatan antara data mentah dan hasil analisis ilmiah.

# Apa itu read, mate, dan coverage?

Sebelum terjun lebih jauh, mari kita luruskan beberapa istilah kunci:

-   ***Read*** - adalah potongan sekuens DNA hasil dari satu kali pembacaan oleh mesin sekuenser.

-   ***Mate*** - mengacu pada pasangan *read* dalam *paired-end* sequencing; satu dari ujung 5', satunya lagi dari 3'.

-   ***Fragment*** - adalah molekul DNA asli yang disekuens.

-   ***Coverage*** - adalah ukuran seberapa dalam sebuah lokasi pada genom dibaca ulang oleh *read*.

**Rumus *coverage*:**

$$
\text{Coverage} = \frac{L \times N}{G}
$$

