---
title: "Quality control hasil *sequencing*"
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

>   **Catatan**: Dalam tutorial ini, kita akan bekerja menggunakan pemrograman Bash dan R dalam lingkungan Linux. Anda dapat menggunakan VS code baik untuk membuat script R maupun menjalankan script Bash tersebut via terminal.

# Daftar isi

-   [Pendahuluan]()
-   [Apa itu read, mate, dan coverage?](#apa-itu-read-mate-dan-coverage)
-   [FASTA & FASTQ: Format dasar NGS dan kualitasnya](#fasta--fastq-format-dasar-ngs-dan-kualitasnya)
-   [SAM & BAM: Dari read ke posisi di genom](#sam--bam-dari-read-ke-posisi-di-genom)
-   [GFF/GTF & BED: Menyematkan makna pada genom](#gffgtf--bed-menyematkan-makna-pada-genom)
-   [VCF: Menyimpan informasi variasi genetik](#vcf-menyimpan-informasi-variasi-genetik)
-   [SRA: Gudangnya data NGS dunia](#sra-gudangnya-data-ngs-dunia)
-   [Memahami format = Menguasai analisis](#memahami-format--menguasai-analisis)

# Pendahuluan