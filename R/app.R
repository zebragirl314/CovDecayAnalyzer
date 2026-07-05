#' Run OriLyt Coverage Decay Analyzer Shiny application. Optionally, load in your parameters along with the app.
#'
#' @param directory, Used to pre-load the directory where your .tab files are stored
#' @param metadata, Used to pre-load sample metadata (must be formatted as csv with columns "sample" and "group")
#'
#' @returns shiny.appobj
#' @export
#'
covApp <- function(directory = "/input", metadata = NULL){
  ##---------- UI ----------
  ui <- page_navbar(
    title = tags$span(
      tags$img(
        src = "https://www.r-project.org/logo/Rlogo.png",
        height = "24px",
        style = "margin-right:8px; vertical-align:middle;"
      ),
      "OriLyt Coverage Decay Analyzer"
    ),
    theme = bs_theme(
      version = 5,
      bootswatch = "flatly",
      primary = "#2C7BB6",
      secondary = "#5A6A7A",
      base_font = font_google("IBM Plex Sans"),
      heading_font = font_google("IBM Plex Mono")
    ),
    navbar_options = navbar_options(bg = "#1a2330"),

    ##--- Tab 1: Configuration & Data Loading ---
    nav_panel(
      title = "Setup",
      icon = icon("sliders"),
      layout_sidebar(
        sidebar = sidebar(
          width = 320,
          bg = "#f0f4f8",

          h5("Input Directory", class = "mt-2 mb-1 fw-bold"),
          textInput(
            "base_dir",
            "Base directory path:",
            value = directory,
            placeholder = "/path/to/data"
          ),

          hr(),
          h5("Genome Parameters", class = "mt-2 mb-1 fw-bold"),
          numericInput(
            "genome_length",
            "Genome length (bp):",
            value = 238080,
            min = 1000
          ),
          numericInput(
            "bin_size",
            "Bin size (bp):",
            value = 1000,
            min = 100,
            max = 10000,
            step = 100
          ),
          numericInput(
            "oriLyt_pos",
            "OriLyt position (bp):",
            value = 150000,
            min = 0
          ),

          hr(),
          h5("Sample Metadata", class = "mt-2 mb-1 fw-bold"),
          p("Paste sample metadata as CSV (columns: sample, group):", class = "text-muted small"),
          tags$textarea(
            id = "sample_metadata_text",
            class = "form-control",
            rows = "14",
            style = "font-family: monospace; font-size: 12px;", ifelse(is.null(metadata), "sample,group
13_FVnoGCV_1_S27,FV-GCV
14_FVnoGCV_2_S28,FV-GCV
15_FVnoGCV_3_S29,FV-GCV
16_FVnoGCV_4_S30,FV-GCV
17_FVnoGCV_5_S31,FV-GCV
18_FVyesGCV_1_S32,FV+GCV
19_FVyesGCV_2_S33,FV+GCV
20_FVyesGCV_3_S34,FV+GCV
21_FVyesGCV_4_S35,FV+GCV
22_FVyesGCV_5_S36,FV+GCV)",
            metadata
          )),

          hr(),
          actionButton("load_data", "Load & Process Data", class = "btn-primary w-100 fw-bold"),
          br(),
          br(),
          uiOutput("load_status")
        ),

        # Main panel: show sample table
        card(
          card_header("Sample Metadata Preview"),
          tableOutput("sample_table")
        ),
        card(
          card_header("File Path Preview"),
          verbatimTextOutput("filepath_preview")
        )
      )
    ),

    ##--- Tab 2: Coverage Plot ---
    nav_panel(
      title = "Coverage Plot",
      icon = icon("chart-line"),
      layout_sidebar(
        sidebar = sidebar(
          width = 280,
          bg = "#f0f4f8",

          h5("Plot Controls", class = "fw-bold mt-2 mb-1"),

          sliderInput(
            "dist_xlim",
            "X-axis limit (bp from OriLyt):",
            min = 5000,
            max = 120000,
            value = 50000,
            step = 1000
          ),

          numericInput(
            "loess_span",
            "LOESS span:",
            value = 0.3,
            min = 0.05,
            max = 1.0,
            step = 0.05
          ),

          checkboxInput("show_points", "Show individual data points", value = TRUE),

          hr(),
          h5("Region Filter (optional)", class = "fw-bold mb-1"),
          checkboxInput("use_region_filter", "Filter to genomic region", value = FALSE),
          conditionalPanel(
            condition = "input.use_region_filter == true",
            numericInput("region_start", "Region start (bp):", value = 135000),
            numericInput("region_end", "Region end (bp):", value = 160000)
          ),

          hr(),
          downloadButton("download_cov_plot", "Download Plot", class = "btn-outline-primary w-100"),
          br(),
          br(),
          downloadButton("download_cov_table", "Download Table (CSV)", class = "btn-outline-secondary w-100")
        ),

        card(
          full_screen = TRUE,
          card_header("Normalized Coverage vs Distance from OriLyt"),
          plotOutput("coverage_plot", height = "520px")
        )
      )
    ),

    ##--- Tab 3: Decay Slopes ---
    nav_panel(
      title = "Decay Slopes",
      icon = icon("arrow-trend-down"),
      layout_sidebar(
        sidebar = sidebar(
          width = 280,
          bg = "#f0f4f8",

          h5("Plot Controls", class = "fw-bold mt-2 mb-1"),
          checkboxInput("show_jitter", "Show individual sample points", value = TRUE),
          checkboxInput("show_outlier", "Show boxplot outliers", value = FALSE),

          hr(),
          h5("Statistical Tests", class = "fw-bold mb-1"),
          selectInput(
            "stat_test",
            "Test to run:",
            choices = c(
              "One-way ANOVA" = "anova",
              "Kruskal-Wallis" = "kruskal",
              "Pairwise Permutation" = "permutation"
            ),
            selected = "anova"
          ),
          actionButton("run_test", "Run Test", class = "btn-primary w-100"),

          hr(),
          downloadButton("download_slope_plot", "Download Plot", class = "btn-outline-primary w-100"),
          br(),
          br(),
          downloadButton("download_slope_table", "Download Slopes (CSV)", class = "btn-outline-secondary w-100")
        ),

        layout_columns(
          col_widths = c(7, 5),
          card(
            full_screen = TRUE,
            card_header("Distribution of Coverage Decay Slopes by Group"),
            plotOutput("slopes_plot", height = "480px")
          ),
          layout_columns(
            col_widths = 12,
            card(
              card_header("Slope Summary by Group"),
              tableOutput("slope_summary_table")
            ),
            card(
              card_header("Statistical Test Results"),
              verbatimTextOutput("stat_test_output")
            )
          )
        )
      )
    ),

    ##--- Tab 4: Per-Sample Explorer ---
    nav_panel(
      title = "Sample Explorer",
      icon = icon("microscope"),
      layout_sidebar(
        sidebar = sidebar(
          width = 260,
          bg = "#f0f4f8",
          h5("Select Sample", class = "fw-bold mt-2 mb-1"),
          uiOutput("sample_selector"),
          hr(),
          h5("Display Options", class = "fw-bold mb-1"),
          checkboxInput("log_scale_sample", "Log-scale Y axis", value = FALSE),
          numericInput(
            "sample_xlim",
            "X-axis max (bp):",
            value = 50000,
            min = 5000
          )
        ),
        card(
          full_screen = TRUE,
          card_header("Per-Sample Coverage Profile"),
          plotOutput("sample_plot", height = "480px")
        )
      )
    )
  )

  ##---------- SERVER ----------

  server <- function(input, output, session) {
    # Reactive: parse sample metadata
    sample_info <- reactive({
      req(input$sample_metadata_text)
      tryCatch(
        read.csv(
          text = input$sample_metadata_text,
          stringsAsFactors = FALSE
        ),
        error = function(e)
          NULL
      )
    })

    output$sample_table <- renderTable({
      sample_info()
    }, striped = TRUE, hover = TRUE, bordered = TRUE)

    output$filepath_preview <- renderText({
      si <- sample_info()
      req(si)
      example_sample <- si$sample[1]
      paste0(
        "Example file path:\n",
        file.path(
          input$base_dir,
          example_sample,
          "08_mutation_identification",
          "Exported.coverage.tab"
        )
      )
    })

    # ---- Data Loading ----
    all_data <- reactiveVal(NULL)
    load_log  <- reactiveVal("")

    observeEvent(input$load_data, {
      si <- sample_info()
      req(si)

      base_dir      <- isolate(input$base_dir)
      bin_size      <- isolate(input$bin_size)
      genome_length <- isolate(input$genome_length)
      oriLyt_pos    <- isolate(input$oriLyt_pos)

      load_log("Loading data...")

      withProgress(message = "Loading samples...", value = 0, {
        n <- nrow(si)
        result_list <- list()
        errors <- c()

        for (i in seq_len(n)) {
          sname <- si$sample[i]
          incProgress(1 / n, detail = paste("Processing:", sname))

          file_path <- file.path(base_dir,
                                 sname,
                                 "08_mutation_identification",
                                 "Exported.coverage.tab")

          tryCatch({
            df <- fread(file_path)
            df[, total_cov := unique_top_cov + unique_bot_cov + redundant_top_cov + redundant_bot_cov]
            df[, bin := ((position - 1) %/% bin_size) + 1]

            binned <- df[, .(mean_cov = mean(total_cov)), by = bin]
            binned[, bin_start := (bin - 1) * bin_size + 1]
            binned[, bin_end   := pmin(bin_start + bin_size - 1, genome_length)]
            binned[, bin_mid   := (bin_start + bin_end) / 2]

            # Normalize to median
            med <- median(binned$mean_cov)
            binned[, norm_cov := mean_cov / med]

            # Distance from OriLyt
            binned[, dist_from_ori := abs(bin_mid - oriLyt_pos)]
            binned[, sample := sname]

            result_list[[i]] <- binned
          }, error = function(e) {
            errors <<- c(errors, paste(sname, ":", conditionMessage(e)))
          })
        }

        combined <- rbindlist(result_list)
        combined  <- merge(combined, si, by = "sample")

        # Log-transform
        epsilon <- 1e-5
        combined[, log_norm_cov := log(norm_cov + epsilon)]

        all_data(combined)

        if (length(errors) > 0) {
          load_log(paste0(
            "Loaded with errors:\n",
            paste(errors, collapse = "\n")
          ))
        } else {
          load_log(paste0("Successfully loaded ", n, " samples."))
        }
      })
    })

    output$load_status <- renderUI({
      msg <- load_log()
      if (msg == "")
        return(NULL)
      cls <- if (grepl("Success", msg))
        "alert alert-success"
      else if (grepl("errors", msg))
        "alert alert-warning"
      else
        "alert alert-info"
      div(class = cls, style = "font-size:13px; white-space:pre-wrap;", msg)
    })

    # ---- Slopes computation ----
    slopes_data <- reactive({
      dat <- all_data()
      req(dat)

      fit <- dat[, broom::tidy(lm(log_norm_cov ~ dist_from_ori)), by = sample]
      slopes <- fit[term == "dist_from_ori", .(sample, slope = estimate)]
      slopes <- merge(slopes, unique(dat[, .(sample, group)]), by = "sample")
      slopes$group <- factor(slopes$group, levels = sort(unique(slopes$group)))
      slopes
    })

    # ---- Coverage Plot ----
    coverage_plot_obj <- reactive({
      dat <- all_data()
      req(dat)

      if (input$use_region_filter) {
        dat <- dat[bin_mid >= input$region_start &
                     bin_mid <= input$region_end]
      }

      dat$group <- factor(dat$group, levels = sort(unique(dat$group)))

      p <- ggplot(dat, aes(x = dist_from_ori, y = norm_cov, color = group))

      if (input$show_points) {
        p <- p + geom_point(alpha = 0.15, size = 0.5)
      }

      p <- p +
        geom_smooth(
          method = "loess",
          span = input$loess_span,
          se = TRUE,
          linewidth = 1.2
        ) +
        coord_cartesian(xlim = c(0, input$dist_xlim)) +
        scale_color_brewer(palette = "Set1") +
        labs(
          title = "Normalized Coverage vs Distance from OriLyt",
          x = "Distance from OriLyt (bp)",
          y = "Normalized Coverage (median-scaled)",
          color = "Group"
        ) +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(face = "bold"),
              legend.position = "top")
      p
    })

    output$coverage_plot <- renderPlot({
      coverage_plot_obj()
    })

    output$download_cov_plot <- downloadHandler(
      filename = "coverage_vs_distance.tiff",
      content  = function(file)
        ggsave(
          file,
          coverage_plot_obj(),
          width = 10,
          height = 6,
          dpi = 300
        )
    )

    output$download_cov_table <- downloadHandler(
      filename = "all_samples_binned_norm_cov.csv",
      content  = function(file)
        fwrite(all_data(), file)
    )

    # ---- Slopes Plot ----
    slopes_plot_obj <- reactive({
      sl <- slopes_data()
      req(sl)

      p <- ggplot(sl, aes(x = group, y = slope, fill = group)) +
        geom_boxplot(alpha = 0.65, outlier.shape = if (input$show_outlier)
          19
          else
            NA) +
        scale_fill_brewer(palette = "Set1") +
        labs(title = "Distribution of Coverage Decay Slopes", x = "Group", y = "Slope (log-linear decay)") +
        theme_minimal(base_size = 14) +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold"))

      if (input$show_jitter) {
        p <- p + geom_jitter(width = 0.2,
                             size = 2.5,
                             alpha = 0.8)
      }
      p
    })

    output$slopes_plot <- renderPlot({
      slopes_plot_obj()
    })

    output$slope_summary_table <- renderTable({
      sl <- slopes_data()
      req(sl)
      sl[, .(
        mean_slope = round(mean(slope), 8),
        sd_slope   = round(sd(slope), 8),
        n          = .N
      ), by = group]
    }, striped = TRUE, hover = TRUE, bordered = TRUE)

    output$stat_test_output <- renderPrint({
      req(input$run_test)
      sl <- slopes_data()
      req(sl)

      isolate({
        if (input$stat_test == "anova") {
          res <- aov(slope ~ group, data = sl)
          print(summary(res))
          cat("\nShapiro-Wilk normality test per group:\n")
          print(by(sl$slope, sl$group, shapiro.test))

        } else if (input$stat_test == "kruskal") {
          res <- kruskal.test(slope ~ group, data = sl)
          print(res)

        } else if (input$stat_test == "permutation") {
          res <- rcompanion::pairwisePermutationTest(slope ~ group,
                                         data = as.data.frame(sl),
                                         distribution = "approximate")
          print(res)
        }
      })
    })

    output$download_slope_plot <- downloadHandler(
      filename = "decay_slopes_boxplot.tiff",
      content  = function(file)
        ggsave(
          file,
          slopes_plot_obj(),
          width = 6,
          height = 8,
          dpi = 300
        )
    )

    output$download_slope_table <- downloadHandler(
      filename = "coverage_decay_slopes.csv",
      content  = function(file)
        fwrite(slopes_data(), file)
    )

    # ---- Sample Explorer ----
    output$sample_selector <- renderUI({
      dat <- all_data()
      req(dat)
      selectInput("selected_sample", "Sample:", choices = unique(dat$sample))
    })

    output$sample_plot <- renderPlot({
      dat <- all_data()
      req(dat, input$selected_sample)

      sdat <- dat[sample == input$selected_sample]
      grp  <- unique(sdat$group)

      y_var <- if (input$log_scale_sample)
        "log_norm_cov"
      else
        "norm_cov"
      y_lab <- if (input$log_scale_sample)
        "Log Normalized Coverage"
      else
        "Normalized Coverage"

      ggplot(sdat, aes_string(x = "dist_from_ori", y = y_var)) +
        geom_line(color = "#2C7BB6", linewidth = 0.8) +
        geom_smooth(
          method = "lm",
          color = "#D7191C",
          linetype = "dashed",
          se = FALSE
        ) +
        coord_cartesian(xlim = c(0, input$sample_xlim)) +
        labs(
          title    = paste0("Coverage Profile: ", input$selected_sample, " (", grp, ")"),
          subtitle = "Red dashed line = log-linear fit",
          x = "Distance from OriLyt (bp)",
          y = y_lab
        ) +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(face = "bold"))
    })
  }

  ##---------- RUN ----------
  shinyApp(ui, server)
}
