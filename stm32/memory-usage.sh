#!/bin/sh

grep '^Sections:' Release/svd_stm32.out.dump -A 18 | grep -v '^        '
