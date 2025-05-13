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

R ada bahasa pemrograman yang didesain khusus untuk analisis data. Untuk menginstalnya, kita bisa merujuk pada halaman instalasi dari [CRAN](https://cloud.r-project.org/). Di Ubuntu, kita bisa mengikuti langkah-langkah berikut.

1.  Buka terminal Ubuntu anda (wsl dengan cara: `wsl -d Ubuntu`), pastikan anda berada di direktori `/home` dan lakukan update sistem dengan menjalankan perintah:

    ```bash
    sudo apt update && sudo apt upgrade
    ```

2.  Install paket pembatu:

    ```bash
    sudo apt install --no-install-recommends software-properties-common dirmngr
    ```

3.  Download repositori resmi dari CRAN untuk memastikan bahwa kita akan terus mendapatkan update terbaru

    ```bash
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    ```

    Perintah ini akan mendownload *public key repository* CRAN untuk Ubuntu (khusus dari maintainer Michael Rutter) menggunakan `wget`, lalu langsung meneruskannya (`|` - *pipe*) ke perintah `tee` untuk disimpan sebagai file kunci baru di sistem, yaitu: `/etc/apt/trusted.gpg.d/cran_ubuntu_key.asc`

    Langkah ini penting supaya sistem mempercayai paket-paket yang diunduh dari CRAN, sehingga proses instalasi software R dari repository tersebut tidak diblokir karena alasan keamanan.

4.  Tambahkan repositori dari CRAN ke `/etc/apt` dengan perintah berikut:

    ```bash
    sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
    ```

    Ini menggunakan akses *superuser* (`sudo`) untuk menambahkan repositori "https://cloud.r-project.org/bin/linux/ubuntu" dan secara otomatis menggunakan versi rilis "*noble*" karena kita menggunakan Ubuntu versi 24. 

5.  Terakhir, install R dengan perintah `sudo`:

    ```bash
    sudo apt install --no-install-recommends r-base
    ```

    Hingga tahap ini, anda sudah berhasil install R di dalam Linux. Untuk cek instalasinya cukup ketikkan `R` pada terminal, maka akan muncul detail versi R yang digunakan.

<div style="position: relative; padding-bottom: 75%; height: 0; overflow: hidden;">
  <iframe src="https://www.youtube.com/embed/f8TM83h4DbA"
          style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"
          frameborder="0"
          allowfullscreen></iframe>
</div>

# Instalasi Anaconda dan tools pendukung

# Instalasi IDE (Integrated development environment)

# Instalasi git