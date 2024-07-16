#!/bin/sh
#本脚本可以将CP2K计算的cp2k-pos-1.xyz中最后一步结构
#转化为cif格式且保留了晶格参数，可直接导入MS和VESTA中
#PS：需要笛卡尔坐标的POSCAR文件


#https://cloud.tencent.com/developer/article/1629932
while getopts ":h" optname
do
    case "$optname" in
      "h")
        echo "#本脚本可以将CP2K计算的cp2k-pos-1.xyz中最后一步结构
#转化为cif格式且保留了晶格参数，可直接导入MS和VESTA中
#PS：需要笛卡尔坐标的POSCAR文件"
	exit
        ;;
      ":")
        echo "No argument value for option $OPTARG"
        ;;
      "?")
        echo "Unknown option $OPTARG"
        ;;
      *)
        echo "Unknown error while processing options"
        ;;
    esac
    #echo "option index is $OPTIND"
done

n0=`sed -n 1p cp2k-pos-1.xyz  | awk '{print $1}'`
n=$(($n0 + 2))
tail -$n cp2k-pos-1.xyz > last.xyz #提取cp2k-pos-1.xyz中最后一步结构

#提取POSCAR的前7行，与最后一步结构的坐标合成新的POSCAR文件last
sed -n '1,7p' POSCAR > last
echo "Cartesian" >> last
sed -n '3,$p' last.xyz | awk '{printf "%16s %16s %16s\n",$2,$3,$4}'>> last

#用vaspkit将last转化为cif
(echo 406;echo 3;echo last; echo 4063) | vaspkit | grep 123

