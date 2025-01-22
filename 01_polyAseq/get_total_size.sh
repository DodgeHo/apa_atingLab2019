#!/bin/bash
# SRA 文件列表
srrs=(SRR4091084 SRR4091085 SRR4091086 SRR4091087 SRR4091088 SRR4091089 SRR4091090 SRR4091091 SRR4091104 SRR4091105 SRR4091106 SRR4091107 SRR4091108 SRR4091109 SRR4091110 SRR4091111 SRR4091113 SRR4091115 SRR4091117 SRR4091119)

# 初始化总大小
total_size=0

# 遍历每个 SRA 文件
for srr in "${srrs[@]}"
do
  # 使用 prefetch 获取文件大小
  size=$(prefetch -s $srr 2>&1 | grep -oP 'size: \K[0-9.]+[GM]')

  # 将大小转换为字节
  if [[ $size == *G ]]; then
    size_bytes=$(echo "${size%G} * 1024 * 1024 * 1024" | bc)
  elif [[ $size == *M ]]; then
    size_bytes=$(echo "${size%M} * 1024 * 1024" | bc)
  else
    size_bytes=$size
  fi

  # 累加总大小
  total_size=$(echo "$total_size + $size_bytes" | bc)

  # 打印当前文件大小
  echo "$srr: $size"
done

# 打印总大小
total_size_gb=$(echo "scale=2; $total_size / 1024 / 1024 / 1024" | bc)
echo "Total size: ${total_size_gb}G"
