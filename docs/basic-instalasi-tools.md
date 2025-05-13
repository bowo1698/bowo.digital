---
title: "Instalasi Tools Bioinformatika"
layout: default
parent: 1. Materi dasar
nav_order: 3
---

# Daftar isi

# Pendahuluan 

Sebelum mempelajari bioinformatika, kita perlu memahami terlebih dahulu apa saja perangkat dan lingkungan kerja yang dibutuhkan untuk menunjang berbagai analisis yang akan dilakukan. Seperti halnya memancing ikan, keberhasilan sangat bergantung pada persiapan dan pemilihan alat yang tepat, dari joran, umpan, hingga teknik memancing. Dalam bioinformatika, "alat pancing" kita adalah software dan lingkungan komputasi, yang masing-masing memiliki fungsi dan keunggulan tersendiri. 

Misalnya dalam bahasa pemrograman, R sangat ungul dalam analisis statistik, karena menawarkan kesederhanaan dan keberlimpahan paket mulai dari yang disediakan oleh `r-base` hingga `Bioconductor`. Sementara itu, Python unggul dalam fleksibilitas dan integrasi dengan pipeline otomatis, serta didukung oleh pustaka seperti `Bioconda`, `Biopython`, `scikit-learn`, dan `pandas`. Adapun bash scripting sangat penting dalam otomasi proses dan pemrosesan file besar di sistem Linux. 

Dalam hal penyejajaran sekuen RNA (RNA-seq alignment) misalnya, tool seperti `STAR` digunakan untuk pemetaan akurat genom referensi. Tool ini sangat cocok untuk studi yang membutuhkan informasi detail tentang struktur gen, terutama dalam mendeteksi titik sambung antar-ekson (*splice junctions*) serta mengidentifikasi varian transkrip mRNA yang dihasilkan melalui *alternative splicing*. Di sisi lain, `Salmon` dan `Kallisto` menggunakan pendekatan *pseudoalignment* yang jauh lebih cepat dan efisien dalam menghitung ekspresi gen, namun tidak mampu mendapatkan informasi tentang struktur gen seperti halnya `STAR`.

Oleh karena itu, sebelum masuk ke tahap analisis, kita perlu terlebih dahulu menyiapkan lingkungan kerja, mulai dari instalasi R, Anaconda, hingga IDE (*integrated development environemnt*) seperti RStudio atau VS Code, yang akan digunakan sebagai "kotak alat" dalam menjalankan berbagai tahap bioinformatika selain menggunakan Terminal bash.

# Instalasi R dan tools pendukung


# Instalasi Anaconda dan tools pendukung

# Instalasi IDE (Integrated development environment)

# Instalasi git