#用excel导入数据, 格式为csv
ori.data <- read.csv("lesson8.csv", header = F)
#以矩阵的方式读入数据, 按行排列, 每三列换一行
data <- matrix(as.matrix(ori.data), nrow(ori.data) / 3, 3, byrow = TRUE)
#关闭区域特定的时间编码方式
Sys.setlocale("LC_TIME", "C")
#用as.POSIXlt()读入字符串数据并转化为date数据, 赋值给date, 或as.Date()
date <- as.POSIXlt(data[, 1], tz = "", "%a %b %d %H:%M:%S HKT %Y")
#对ip和pv所在的列转化为数值型
IP <- as.numeric(data[, 2])
PV <- as.numeric(data[, 3])
head(data)
#恢复区域特地的时间编码方式
Sys.setlocale("LC_TIME", "")
#用ggplot2绘图
library(ggplot2)
#用reshape包中的melt函数分解数据
library(reshape2)
p.data <- data.frame(date, IP, PV)
meltdata <- melt(p.data, id = (c("date")))
#用对IP和PV做分页处理, y轴刻度自由变化
graphic <- ggplot(data = meltdata, aes(x = date, y = value, color = variable)) + geom_line() + geom_point()
graphic <- graphic + facet_grid(variable ~ ., scales = "free_y")
#美化, 添加标题, 坐标, 更改图例
graphic<- graphic + labs(x = "日期", y = "人次", title = "某网站7月至10月IP/PV统计") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  scale_colour_discrete(name = "",labels = c("IP","PV")) +
  theme(strip.text.y = element_text(angle = 0))