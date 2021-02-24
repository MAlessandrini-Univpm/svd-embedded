# svd-embedded
SVD algorithms for embedded systems (STM32 MCUs, Cortex-M).

This is the software used in the article:

Alessandrini, M.; Biagetti, G.; Crippa, P.; Falaschetti, L.; Manoni, L.; Turchetti, C. "Singular Value Decomposition in Embedded Systems Based on ARM Cortex-M Architecture". Electronics 2021, 10, 34. https://doi.org/10.3390/electronics10010034 

If you use this work, please cite the article:
```
@Article{electronics10010034,
AUTHOR = {Alessandrini, Michele and Biagetti, Giorgio and Crippa, Paolo and Falaschetti, Laura and Manoni, Lorenzo and Turchetti, Claudio},
TITLE = {Singular Value Decomposition in Embedded Systems Based on ARM Cortex-M Architecture},
JOURNAL = {Electronics},
VOLUME = {10},
YEAR = {2021},
NUMBER = {1},
ARTICLE-NUMBER = {34},
URL = {https://www.mdpi.com/2079-9292/10/1/34},
ISSN = {2079-9292},
ABSTRACT = {Singular value decomposition (SVD) is a central mathematical tool for several emerging applications in embedded systems, such as multiple-input multiple-output (MIMO) systems, data analytics, sparse representation of signals. Since SVD algorithms reduce to solve an eigenvalue problem, that is computationally expensive, both specific hardware solutions and parallel implementations have been proposed to overcome this bottleneck. However, as those solutions require additional hardware resources that are not in general available in embedded systems, optimized algorithms are demanded in this context. The aim of this paper is to present an efficient implementation of the SVD algorithm on ARM Cortex-M. To this end, we proceed to (i) present a comprehensive treatment of the most common algorithms for SVD, providing a fairly complete and deep overview of these algorithms, with a common notation, (ii) implement them on an ARM Cortex-M4F microcontroller, in order to develop a library suitable for embedded systems without an operating system, (iii) find, through a comparative study of the proposed SVD algorithms, the best implementation suitable for a low-resource bare-metal embedded system, (iv) show a practical application to Kalman filtering of an inertial measurement unit (IMU), as an example of how SVD can improve the accuracy of existing algorithms and of its usefulness on a such low-resources system. All these contributions can be used as guidelines for embedded system designers. Regarding the second point, the chosen algorithms have been implemented on ARM Cortex-M4F microcontrollers with very limited hardware resources with respect to more advanced CPUs. Several experiments have been conducted to select which algorithms guarantee the best performance in terms of speed, accuracy and energy consumption.},
DOI = {10.3390/electronics10010034}
}
```
# Usage

The code has been tested on a 32F429IDISCOVERY board, but can be used with other STM32 devices, provided the memory requirements are met. See the article for details on the hardware setup.

The code is designed to be compiled in command-line on a GNU/Linux system. From a shell in the stm32 directory, type `make` to compile the code. You must have the `gcc-arm-none-eabi` suite installed.

To flash and interactively debug the code, you can use the `make dbg` command, that uses gdb and openocd to communicate with the board.
