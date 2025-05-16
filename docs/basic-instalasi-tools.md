---
title: "Instalasi Tools Bioinformatika"
layout: default
parent: 1. Materi dasar
nav_order: 5
---

<p style="text-align: right; font-size: 0.9rem;">
  <a href="https://www.bowo.digital/" style="font-weight: bold;">← Beranda</a>
  &nbsp;&nbsp;|&nbsp;&nbsp;
  <a href="https://www.bowo.digital/docs/part1.html" style="font-weight: bold;">Pilih materi →</a>
</p>

<h1 style="text-align: center; font-size: 2.5rem; font-weight: bold; margin-bottom: 0.5rem;">
  <a href="https://www.bowo.digital/docs/basic-instalasi-tools.html" style="text-decoration: none; color: inherit;">
    Instalasi Tools Bioinformatika
  </a>
</h1>

<p style="text-align: center; font-size: 1.2rem;">
  Oleh <a href="https://www.bowo.digital/docs/bio.html" target="_blank">Agus Wibowo</a>
</p>

<div style="text-align: center; margin-bottom: 1.5rem;">
  <img src="https://static.vecteezy.com/system/resources/thumbnails/055/370/498/small_2x/scientist-analyzing-dna-strands-on-a-digital-screen-with-glowing-data-in-a-lab-environment-photo.jpg" alt="Tools instalation" style="max-width: 100%; height: auto;">
</div>

# Daftar isi

-   [Pendahuluan](#pendahuluan)
-   [Instalasi R](#instalasi-r)
-   [Instalasi Anaconda](#instalasi-anaconda)
-   [Instalasi IDE](#instalasi-ide-integrated-development-environment)
-   [Basic instalasi tools R dan Conda](#basic-instalasi-tools-r-dan-conda)

# Pendahuluan 

Sebelum mempelajari bioinformatika, kita perlu memahami terlebih dahulu apa saja perangkat dan lingkungan kerja yang dibutuhkan untuk menunjang berbagai analisis yang akan dilakukan. Seperti halnya memancing ikan, keberhasilan sangat bergantung pada persiapan dan pemilihan alat yang tepat, dari joran, umpan, hingga teknik memancing. Dalam bioinformatika, "alat pancing" kita adalah software dan lingkungan komputasi, yang masing-masing memiliki fungsi dan keunggulan tersendiri. 

Misalnya dalam bahasa pemrograman, R sangat ungul dalam analisis statistik, karena menawarkan kesederhanaan dan keberlimpahan paket mulai dari yang disediakan oleh `r-base` hingga `Bioconductor`. Sementara itu, Python unggul dalam fleksibilitas dan integrasi dengan pipeline otomatis, serta didukung oleh pustaka seperti `Bioconda`, `Biopython`, `scikit-learn`, dan `pandas`. Adapun bash scripting sangat penting dalam otomasi proses dan pemrosesan file besar di sistem Linux. 

Dalam hal penyejajaran sekuen RNA (RNA-seq alignment) misalnya, tool seperti `STAR` digunakan untuk pemetaan akurat genom referensi. Tool ini sangat cocok untuk studi yang membutuhkan informasi detail tentang struktur gen, terutama dalam mendeteksi titik sambung antar-ekson (*splice junctions*) serta mengidentifikasi varian transkrip mRNA yang dihasilkan melalui *alternative splicing*. Di sisi lain, `Salmon` dan `Kallisto` menggunakan pendekatan *pseudoalignment* yang jauh lebih cepat dan efisien dalam menghitung ekspresi gen, namun tidak mampu mendapatkan informasi tentang struktur gen seperti halnya `STAR`.

Oleh karena itu, sebelum masuk ke tahap analisis, kita perlu terlebih dahulu menyiapkan lingkungan kerja, mulai dari instalasi R, Anaconda, hingga IDE (*integrated development environemnt*) seperti RStudio atau VS Code, yang akan digunakan sebagai "kotak alat" dalam menjalankan berbagai tahap bioinformatika selain menggunakan Terminal bash.

# Instalasi R

R ada bahasa pemrograman yang didesain khusus untuk analisis data. Untuk menginstalnya, kita bisa merujuk pada halaman instalasi dari [CRAN](https://cloud.r-project.org/). Di Ubuntu atau UNIX lainnya (kali ini saya menggunakan Ubuntu via WSL), kita bisa mengikuti langkah-langkah berikut.

1.  Buka terminal Ubuntu anda (wsl dengan cara: `wsl -d Ubuntu`), pastikan anda berada di direktori `/home` dan lakukan update sistem dengan menjalankan perintah:

    ```bash
    # bash
    sudo apt update && sudo apt upgrade
    ```

2.  Install paket pembatu:

    ```bash
    # bash
    sudo apt install --no-install-recommends software-properties-common dirmngr
    ```

3.  Download repositori resmi dari CRAN untuk memastikan bahwa kita akan terus mendapatkan update terbaru

    ```bash
    # bash
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    ```

    Perintah ini akan mendownload *public key repository* CRAN untuk Ubuntu (khusus dari maintainer Michael Rutter) menggunakan `wget`, lalu langsung meneruskannya (`|` - *pipe*) ke perintah `tee` untuk disimpan sebagai file kunci baru di sistem, yaitu: `/etc/apt/trusted.gpg.d/cran_ubuntu_key.asc`

    Langkah ini penting supaya sistem mempercayai paket-paket yang diunduh dari CRAN, sehingga proses instalasi software R dari repository tersebut tidak diblokir karena alasan keamanan.

4.  Tambahkan repositori dari CRAN ke `/etc/apt` dengan perintah berikut:

    ```bash
    # bash
    sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
    ```

    Ini menggunakan akses *superuser* (`sudo`) untuk menambahkan repositori "https://cloud.r-project.org/bin/linux/ubuntu" dan secara otomatis menggunakan versi rilis "*noble*" karena kita menggunakan Ubuntu versi 24. 

5.  Terakhir, install R dengan perintah `sudo`:

    ```bash
    # bash
    sudo apt install --no-install-recommends r-base
    ```

    Hingga tahap ini, anda sudah berhasil install R di dalam Linux. Untuk cek instalasinya cukup ketikkan `R` pada terminal, maka akan muncul interpreter R dan detail versi yang digunakan.

<div style="position: relative; padding-bottom: 75%; height: 0; overflow: hidden;">
  <iframe src="https://www.youtube.com/embed/f8TM83h4DbA"
          style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"
          frameborder="0"
          allowfullscreen></iframe>
</div>

# Instalasi Anaconda

Anaconda (bukan ular) merupakan distribusi *open-source* dari bahasa pemrograman Python yang dirancang khusus untuk *data science*, *machine learning*, dan *scientific computing*. Ditrubusi ini mencakup: Python interpreter, Conda manager, dan ribuan paket populer seperti Numpy, scikit-learn, hingga paket-paket bioinformatika seperti samtools, RSEM, STAR, HISAT, dan masih banyak lagi. Untuk instalasinya di Linux, kita bisa mengikuti tutorial berikut.

1.  Unduh installer Anaconda

    Kunjungi situs resmi Anaconda [di sini](https://www.anaconda.com/download/success) dan salin tautan unduhan untuk Linux dengan cara memilih distribusi `Linux - <arsitektur> Installer` yang sesuai dengan arsitektur komputer anda. Misal karena saya menggunakan ARM64, maka saya pilih ARM64. Kemudian, klik kanan dan pilih "*Copy Link Address*" untuk mendapatkan tautan unduhan. Buka terminal bash dan gunakan `wget`

    ```bash
    # bash
    wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-aarch64.sh
    ```

    File hasil download akan disimpan pada direktori dimana Anda menjalankan Bash, misalnya Anda menjalankan Bash dari direktori `/home/bowo/Downloads`, maka file `Anaconda3-2024.10-1-Linux-aarch64.sh` akan disimpan di direktori ini.

2.  Jalankan installer

    Karena file installer yang kita unduh ber-ekstensi `.bash`, maka ini bisa diinstal dengan menggunakan Bash secara langsung.

    ```bash
    # bash
    bash Anaconda3-2024.10-1-Linux-aarch64.sh
    ```

    Ikuti semua petunjuk yang muncul. Anda hanya perlu tekan tombol `ENTER + SPACE` untuk langsung menuju bagian akhir dimana anda harus ketik "yes" untuk menerima *license terms*. Di akhir, tekan ENTER untuk menginstall `base` PREFIX secara default di terminal, sehingga Anacona dapat dipanggil hanya dengan mengetikkan `conda`. 

<div style="position: relative; padding-bottom: 75%; height: 0; overflow: hidden;">
  <iframe src="https://www.youtube.com/embed/LzuMruPOT9c"
          style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"
          frameborder="0"
          allowfullscreen></iframe>
</div>

# Instalasi IDE (*Integrated development environment*)

Agar dapat menulis, menjalankan, dan mengelola skrip bioinformatika secara efisien, kita memerlukan sebuah IDE (Integrated Development Environment) atau lingkungan pengembangan terpadu. Dalam semua tutorial di situs ini, kita akan menggunakan Visual Studio Code (VS Code) sebagai IDE utama.

VS Code dipilih karena bersifat ringan, fleksibel, dan mudah diintegrasikan dengan berbagai kernel interpreter seperti R, Python, dan Bash, serta mendukung lingkungan Conda secara langsung. Hal ini menjadikan VS Code sebagai pusat kendali yang ideal untuk seluruh proses analisis.

Lebih lanjut, dengan dukungan ekstensi seperti Jupyter, kita dapat menjalankan skrip-skrip tersebut dalam bentuk notebook interaktif, sehingga kode, output, dan dokumentasi dapat ditulis dan dibaca dalam satu tampilan yang terstruktur dan mudah dikelola.

Instalasi VS code dan semua pengaturannya dapat mengikuti panduan berikut:

1.  Download VS code dari situs resminya [di sini](https://code.visualstudio.com/Download#). Sesuikan dengan arsitektur komputer yang anda miliki. Misalnya disini saya menggunakan Debian dengan basis ARM64, maka saya mendownload versi Linux debian (*.deb) dengan arsitektur ARM64.

2.  Buka terminal dan jalankan perintah berikut satu persatu

    ```bash
    # bash
    cd Downloads/
    chmod +x code_1.100.1-1746807040_arm64.deb
    sudo dpkg -i code_1.100.1-1746807040_arm64.deb
    ```

    Perintah di atas pada dasarnya hanya membuka direktori tempat VS code diunduh (`~/Downloads`), kemudian mengubah hak akses installer VS code supaya dapat dieksekusi dengan perintah `chmod +x` dan menginstal VS code dengan hak akses *superuser* menggunakan perintah `sudo dpkg -i <nama_file_VS_code_installer.deb>`

3.  Untuk mengintegrasikan antara VS code dengan R dan Anaconda, kita dapat menjalankan perintah berikut di dalam interpreter R.

    ```r
    # r
    install.packages(c("languageserver", "IRkernel"))
    IRkernel::installspec()
    ```

4.  Buka VS code dan pada tab Extension (Ctrl + Shift + X) install extension diantaranya: R, Python, dan Jupyter

5.  Hingga langkah ini, kita sudah berhasil menginstall VS code dan mengintegrasikannya dengan R dan Anaconda.


# Basic instalasi tools R dan Conda

Untuk instalasi tools, baik R maupun Conda memiliki cara yang berbeda.

Perintah dasar instalasi pada R secara umum untuk kebutuhan bioinformatika, dibagi menjadi 3 jenis:

1.  Instalasi tools yang berasal dari r-base CRAN

    Semua tools yang berasal dari CRAN yang biasanya untuk kebutuhan analisis data statistik, dapat diinstal dengan menggunakan perintah `install.packages("nama_tools")` di dalam R interpreter:

    ```r
    # r
    # misalnya tidyverse
    install.packages("tidyverse")

    # atau jika menginstal lebih dari satu tools bisa menggunakan fungsi c()
    install.packages(c("tools_1", "tools_2", "tools3"))

    # misalnya
    install.packages(c("gert","usethis"))
    ```

2.  Instalasi tools yang berasal dari Bioconductor

    Banyak tools bioinformatika seperti `DESeq2`, `limma`, `GenomicFeatures`, dll, hanya tersedia melalui Bioconductor. Jadi kita perlu menginstall terlebih dahulu `BiocManager` di dalam R, dengan cara sebagai berikut.

    ```r
    # r
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install()
    ```

    Setelah itu, kita dapat menjalankan instalasi tools menggunakan perintah `BiocManager::install("nama_tools")`

    ```r
    # r
    # install DESeq2
    BiocManager::install("DESeq2")

    # Jika tools lebih dari satu, gunakan fungsi c()
    BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "GenomicAlignments", "ComplexHeatmap"))
    ```

3.  Instalasi tools external

    Ada beberapa tools yang tidak masuk ke dalam paket CRAN maupun Bioconductor, biasanya mereka adalah tools yang berasal dari Github dan dapat kita instal dengan perintah `devtools::install_github("nama_tools")`

    ```r
    # r
    # install clusterProfiler dari Github
    devtools::install_github("YuLab-SMU/clusterProfiler")
    ```

**Catatan**

Ketika instalasi paket-paket R, mungkin kita akan menjumpai error atau kegagalan instal lantaran ketiadaan beberapa paket di dalam sistem Linux. Jika demikan, maka kita harus menginstal paket-paket ini satu persatu (lihat video cara mencari paket yang dibutuhkan). Misal dalam kasus saya, `tidyverse` dan `devtools` sempat gagal diinstal karena tidak adanya beberapa paket di Linux saya seperti: `libcurl4-openssl-dev` dan `libgit2-dev`, maka kita harus menginstal paket-paket tersebut sebagai contoh berikut:

```bash
# bash
sudo apt install -y libfreetype6-dev libharfbuzz-dev libfribidi-dev libpng-dev libtiff5-dev libjpeg-dev libssl-dev libfontconfig1-dev libcurl4-openssl-dev libgit2-dev
```

Setelah selesai, ulang lagi instalasi `tidyverse` dan `devtools` dengan menjalankan:

```r
# r
install.packages(c("tidyverse","devtools"))
```

Sementara itu untuk Conda, kita dapat menginstall tools yang dibutuhkan dengan perintah dasar `conda install nama_tools` melalui Bash:

```bash
# bash
# misal menginstal scikit-learn
conda install scikit-learn
```

Namun untuk keperluan bioinformatika, kita perlu mengaktifkan Bioconda, yaitu sebuah *channel* (saluran distribusi paket) yang secara khusus menyediakan ribuan paket bioinformatika siap pakai. Agar Bioconda dapat digunakan, kita harus menambahkan dan mengurutkan *channel*-nya dengan benar menggunakan perintah berikut:

```bash
# bash
# jalankan perintah berikut satu persatu
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Setelah itu, kita bisa membuat beberapa *Environement* di dalam Conda untuk memisahkan setiap tools yang kita instal berdasarkan tujuannya. Misalnya, kita bisa menyatukan semua tools terkait *alignment* didalam satu *environment* bernama `alignemnt_tools`, dan kemudian menginstal tools *alignment* di dalamnya.

```bash
# bash
# buat environment bernama alignment_tools dengan menggunakan Python versi 3.11 
conda create --name alignment_tools python=3.11

# aktifkan environment alignment_tools
conda activate alignment_tools

# instal tools alignment 
conda install samtools bowtie2 hisat2 star bwa salmon
```

Dengan membuat *environment* secara terpisah dengan nama unik, ini bisa membantu kita untuk mengorganisir berbagai proyek atau analisis bioinformatika agar tetap rapi dan terisolasi. Setiap *environment* memiliki versi Python dan dependensi paketnya sendiri, sehingga mencegah konflik antar tools yang berbeda kebutuhan. 

Misalnya, kita bisa membuat *environment*:

-   `rna_seq` untuk analisis ekspresi gen dengan `DESeq2`, `Salmon`, dan `HISAT2`

-   `variant_calling` untuk pipeline analisis varian dengan `BWA`, `GATK`, dan `bcftools`

-   `machine_learning` untuk *genomic prediction* menggunakan `scikit-learn` atau `TensorFlow`

Dengan pendekatan ini, kita tidak hanya menjaga kestabilan sistem, tetapi juga memudahkan reproduksibilitas analisis. Jika ingin melakukan pemindahan, cukup export *environment*-nya menggunakan:

```bash
# bash
# misal ingin mengexport alignment_tools
# aktifkan environment alignment_tools
conda activate alignment_tools

# export
conda env export > alignment_tools.yml
```

<div style="position: relative; padding-bottom: 75%; height: 0; overflow: hidden;">
  <iframe src="https://www.youtube.com/embed/e9EpqCcru5M"
          style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"
          frameborder="0"
          allowfullscreen></iframe>
</div>
```
<p style="text-align: right; font-size: 0.9rem;">
  <a href="https://www.bowo.digital/" style="font-weight: bold;">← Beranda</a>
  &nbsp;&nbsp;|&nbsp;&nbsp;
  <a href="https://www.bowo.digital/docs/part1.html" style="font-weight: bold;">Pilih materi →</a>
</p>