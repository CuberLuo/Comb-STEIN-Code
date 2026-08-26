## 重复试验模拟和汇总聚合。

## 在一个命名场景下重复模拟每种设计方法。
## 返回方法级运行特征表，并将真实 OBDC 和延迟设置保存为属性。
run.scenario <- function(
  name,
  ntrial = DEFAULT_NTRIAL,
  delayed = FALSE,
  delay.control = NULL,
  methods = c(
    "Comb-STEIN-L",
    "Comb-STEIN-H",
    "Runin-STEIN",
    "Waterfall-L",
    "Waterfall-H"
  )
) {
  delayed <- isTRUE(delayed)
  delay.control <- normalize.delay.control(delay.control)

  scenario <- build.pTpE(name)
  pT <- scenario$pT
  pE <- scenario$pE
  true.ij <- get.true.obdc(pT, pE)

  out <- data.frame(
    Scenario = name,
    Method = methods,
    Delay = if (delayed) "TITE" else "Complete",
    Selection = NA,
    Patient = NA,
    Sample = NA,
    DLTs = NA,
    Tox = NA,
    Eff = NA,
    Duration = NA
  )

  for (method in methods) {
    selection.count <- 0
    patient.percent.sum <- 0
    sample.sum <- 0
    dlt.sum <- 0
    toxicity.stop.count <- 0
    efficacy.stop.count <- 0
    duration.sum <- 0

    # 进行独立的ntrial次试验，每次试验运行5种方法
    for (trial in seq_len(ntrial)) {
      result <- tryCatch(
        switch(
          method,
          "Comb-STEIN-L" = simulate.trial.comb.stein(
            pT, pE, "CombSTEIN-L", delayed, delay.control
          ),
          "Comb-STEIN-H" = simulate.trial.comb.stein(
            pT, pE, "CombSTEIN-H", delayed, delay.control
          ),
          "Runin-STEIN" = simulate.trial.runin.stein(
            pT, pE, delayed, delay.control
          ),
          "Waterfall-L" = simulate.trial.waterfall(
            pT, pE, "Waterfall-L", delayed, delay.control
          ),
          "Waterfall-H" = simulate.trial.waterfall(
            pT, pE, "Waterfall-H", delayed, delay.control
          )
        ),
        error = function(err) {
          stop(
            sprintf(
              "Scenario %s, method %s, trial %d failed: %s",
              name,
              method,
              trial,
              conditionMessage(err)
            ),
            call. = FALSE
          )
        }
      )

      ## result 是一次试验的完整输出，主要包含：
      ## n：每个剂量组合实际入组的患者数；
      ## x：每个剂量组合观察到的毒性事件数；
      ## sel：试验结束时选择的剂量组合；
      ## earlystop_tox / earlystop_eff：是否因毒性/疗效原因提前停止；
      ## duration：完整结局或 TITE 模式下的试验持续时间。
      total.n <- sum(result$n)

      ## 累加当前试验的总样本量和总毒性事件数。
      ## 循环结束后，再除以 ntrial 得到平均值。
      sample.sum <- sample.sum + total.n
      dlt.sum <- dlt.sum + sum(result$x)

      ## 统计该次试验是否因毒性或疗效规则提前结束。
      if (result$earlystop_tox) {
        toxicity.stop.count <- toxicity.stop.count + 1
      }
      if (result$earlystop_eff) {
        efficacy.stop.count <- efficacy.stop.count + 1
      }
      duration.sum <- duration.sum + result$duration

      ## 判断最终选中的组合是否属于真实 utility 最大的 OBDC。
      ## true.ij 可能包含多个并列最优组合，因此逐行比较；
      ## result$sel 是当前试验选择的一个 (行, 列) 坐标。
      if (!any(is.na(result$sel))) {
        correct <- any(apply(
          true.ij,
          1,
          function(ij) all(result$sel == ij)
        ))
        if (correct) selection.count <- selection.count + 1
      }

      ## 计算有多少患者被分配到了真实 OBDC。
      ## 如果存在多个并列 OBDC，就把这些组合的患者数全部相加。
      ## 这里使用总入组人数作分母，得到当前试验的百分比。
      if (total.n > 0) {
        true.index <- cbind(true.ij[, 1], true.ij[, 2])
        n.at.true <- sum(result$n[true.index])
        patient.percent.sum <-
          patient.percent.sum + 100 * n.at.true / total.n
      }
    }

    ## 当前方法的所有重复试验完成后，将累计值转换成汇总指标：
    ## Selection：正确选择 OBDC 的试验百分比；
    ## Patient：分配到真实 OBDC 的患者百分比；
    ## Sample：平均入组人数；
    ## DLTs：平均毒性事件数；
    ## Tox / Eff：因毒性/疗效提前停止的试验百分比；
    ## Duration：平均试验持续时间（天）。
    row <- out$Method == method
    out[row, "Selection"] <- round(100 * selection.count / ntrial, 1)
    out[row, "Patient"] <- round(patient.percent.sum / ntrial, 1)
    out[row, "Sample"] <- round(sample.sum / ntrial, 1)
    out[row, "DLTs"] <- round(dlt.sum / ntrial, 1)
    out[row, "Tox"] <- round(100 * toxicity.stop.count / ntrial, 1)
    out[row, "Eff"] <- round(100 * efficacy.stop.count / ntrial, 1)
    out[row, "Duration"] <- round(duration.sum / ntrial, 1)
  }

  ## 将真实 OBDC 和延迟设置作为属性附加到结果表。
  ## 这样 main.R 可以读取这些信息并打印真实 OBDC，
  ## 但它们不会被当作普通数据列写入汇总表。
  attr(out, "true.obdcs") <- true.ij
  attr(out, "delayed") <- delayed
  attr(out, "delay.control") <- delay.control
  out
}
