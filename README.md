102下 陳玉彬熱傳 電腦作業   
===========   
   
題目說明:
-------------  
請閱讀本文件: [computer_program_Q.pdf](https://github.com/janelin612/HeatTransfer/blob/master/computer_program_Q.pdf)   
   
程式碼說明:   
-------------   
本程式支援的網格數量底邊(length)為7的整數倍，高(hight)為3的整數倍(!=3)  
可自行於程式碼的[此處](https://github.com/janelin612/HeatTransfer/blob/master/code/code.m#L29)進行修改    
惟須注意一點為輸入值須自行加一，如7X6則須輸入:
```    
numberOfLength=8;   
numberOfHight=7;
```
   
程式輸出:
-----------
其將會輸出四張圖，分別為網格解析度   
+ 2cmx2cm
+ 2cmx1cm
+ 1cmx2cm
+ 1cmx1cm   
的鰭片溫度分布以及熱通量，如下圖所示:   
![figure](https://raw.githubusercontent.com/janelin612/HeatTransfer/master/figure&data/figure.png)