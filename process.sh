#!usr/bin/bash

input_file=$1
output_dir=$2
username=$3
option=$4
restriction=$5
base_dir=${input_file%/*}

export ARNIEFILE="/var/www/platform/app/scripts/sjqprimer/example_arnie_file.txt"

jobname=${output_dir##*/}
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

cd ${input_file%/*}
set -e

trap "/usr/local/anaconda/bin/python -u /var/www/platform/app/scripts/ngs/taskfailed.py $jobname $username" EXIT
echo "input_file: $input_file" > "$script_dir/args.txt"
echo "output_dir: $output_dir" >> "$script_dir/args.txt"
echo "username: $username" >> "$script_dir/args.txt"
echo "option: $option" >> "$script_dir/args.txt"
echo "restriction: $restriction" >> "$script_dir/args.txt"

echo "参数已写入到 $script_dir/args.txt"


if [ "$option" = "H. sapiens" ];then
	c2w="/var/www/platform/app/scripts/resources/codon/human/compodynamics_human_codon_usage.tsv.usage_per_aa.rare.csv.C2W.csv"
elif [ "$option" = "M. musculus" ];then
	echo "尚无小鼠密码子参考"
	exit 1
elif [ "$option" = "E. coli" ];then
	c2w="/var/www/platform/app/scripts/resources/codon/E_coli/E_coli.compodynamics.codon_usage.tsv.usage_per_aa.tmp.tsv.C2W.csv"
elif [ "$option" = "O. sativa" ];then
	c2w="/var/www/platform/app/scripts/resources/codon/rice/compodynamics_rice_codon_usage.tsv.usage_per_aa.rare.csv.C2W.csv"
elif [ "$option" = "Z. mays" ];then
	c2w="/var/www/platform/app/scripts/resources/codon/maize/compodynamics_maize_codon_usage.tsv.usage_per_aa.rare.csv.C2W.csv"
else
	echo "选项不存在！"
fi

sed '1d' $input_file | awk -F "," -v OFS="\t" '{print ">"$1"\n"$3}' > tmp.fa
finetune_CDS.py -i tmp.fa --C2W $c2w --CAI_s 100 --GC_s 0 --GC3_s 1 --minT_s 5 -r `echo $restriction|awk -F '[, ]+' -v OFS=" " '{ for (i = 1; i <= NF; i++) printf "%s%s", $i, (i==NF ? "\n" : " ") }'`
/var/www/platform/app/scripts/codon_opt/reroll.py $input_file tmp.fa.opt.csv $c2w "`echo $restriction|awk -F '[, ]+' -v OFS=" " '{ for (i = 1; i <= NF; i++) printf "%s%s", $i, (i==NF ? "\n" : " ") }'`" $output_dir
out_parent_dir="/var/www/platform/public/tmp/codon_opt/result"
cd $out_parent_dir
cp ${base_dir}/log.* ${output_dir}/
zip -r "$jobname.zip" $jobname
# 更新数据库
/usr/local/anaconda/bin/python '/var/www/platform/app/scripts/ngs/updateMysql.py' $jobname $username

trap - EXIT
