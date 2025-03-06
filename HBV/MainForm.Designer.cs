namespace HBV
{
    partial class MainForm
    {
        /// <summary>
        ///  Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        ///  Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        ///  Required method for Designer support - do not modify
        ///  the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            components = new System.ComponentModel.Container();
            plotRainfall = new ScottPlot.FormsPlot();
            splitContainer1 = new SplitContainer();
            butRunModel = new Button();
            butRunOptimizer = new Button();
            butLoadData = new Button();
            tabControl1 = new TabControl();
            tabPage1 = new TabPage();
            tabPage2 = new TabPage();
            plotRiverDischarge = new ScottPlot.FormsPlot();
            tmr = new System.Windows.Forms.Timer(components);
            butSave = new Button();
            dlgSave = new SaveFileDialog();
            ((System.ComponentModel.ISupportInitialize)splitContainer1).BeginInit();
            splitContainer1.Panel1.SuspendLayout();
            splitContainer1.Panel2.SuspendLayout();
            splitContainer1.SuspendLayout();
            tabControl1.SuspendLayout();
            tabPage1.SuspendLayout();
            tabPage2.SuspendLayout();
            SuspendLayout();
            // 
            // plotRainfall
            // 
            plotRainfall.Dock = DockStyle.Fill;
            plotRainfall.Location = new Point(3, 3);
            plotRainfall.Margin = new Padding(4, 3, 4, 3);
            plotRainfall.Name = "plotRainfall";
            plotRainfall.Size = new Size(1000, 727);
            plotRainfall.TabIndex = 0;
            // 
            // splitContainer1
            // 
            splitContainer1.Dock = DockStyle.Fill;
            splitContainer1.Location = new Point(0, 0);
            splitContainer1.Name = "splitContainer1";
            // 
            // splitContainer1.Panel1
            // 
            splitContainer1.Panel1.Controls.Add(butSave);
            splitContainer1.Panel1.Controls.Add(butRunModel);
            splitContainer1.Panel1.Controls.Add(butRunOptimizer);
            splitContainer1.Panel1.Controls.Add(butLoadData);
            // 
            // splitContainer1.Panel2
            // 
            splitContainer1.Panel2.Controls.Add(tabControl1);
            splitContainer1.Size = new Size(1184, 761);
            splitContainer1.SplitterDistance = 166;
            splitContainer1.TabIndex = 1;
            // 
            // butRunModel
            // 
            butRunModel.Location = new Point(12, 79);
            butRunModel.Name = "butRunModel";
            butRunModel.Size = new Size(129, 61);
            butRunModel.TabIndex = 2;
            butRunModel.Text = "Run Model";
            butRunModel.UseVisualStyleBackColor = true;
            butRunModel.Click += butRunModel_Click;
            // 
            // butRunOptimizer
            // 
            butRunOptimizer.Location = new Point(12, 146);
            butRunOptimizer.Name = "butRunOptimizer";
            butRunOptimizer.Size = new Size(129, 61);
            butRunOptimizer.TabIndex = 1;
            butRunOptimizer.Text = "Run Optimizer";
            butRunOptimizer.UseVisualStyleBackColor = true;
            butRunOptimizer.Click += butRunOptimizer_Click;
            // 
            // butLoadData
            // 
            butLoadData.Location = new Point(12, 12);
            butLoadData.Name = "butLoadData";
            butLoadData.Size = new Size(129, 61);
            butLoadData.TabIndex = 0;
            butLoadData.Text = "Load and Plot Data";
            butLoadData.UseVisualStyleBackColor = true;
            butLoadData.Click += butLoadAndPlotData_Click;
            // 
            // tabControl1
            // 
            tabControl1.Controls.Add(tabPage1);
            tabControl1.Controls.Add(tabPage2);
            tabControl1.Dock = DockStyle.Fill;
            tabControl1.Location = new Point(0, 0);
            tabControl1.Name = "tabControl1";
            tabControl1.SelectedIndex = 0;
            tabControl1.Size = new Size(1014, 761);
            tabControl1.TabIndex = 1;
            // 
            // tabPage1
            // 
            tabPage1.Controls.Add(plotRainfall);
            tabPage1.Location = new Point(4, 24);
            tabPage1.Name = "tabPage1";
            tabPage1.Padding = new Padding(3);
            tabPage1.Size = new Size(1006, 733);
            tabPage1.TabIndex = 0;
            tabPage1.Text = "rainfall";
            tabPage1.UseVisualStyleBackColor = true;
            // 
            // tabPage2
            // 
            tabPage2.Controls.Add(plotRiverDischarge);
            tabPage2.Location = new Point(4, 24);
            tabPage2.Name = "tabPage2";
            tabPage2.Padding = new Padding(3);
            tabPage2.Size = new Size(1006, 733);
            tabPage2.TabIndex = 1;
            tabPage2.Text = "river discharge";
            tabPage2.UseVisualStyleBackColor = true;
            // 
            // plotRiverDischarge
            // 
            plotRiverDischarge.Dock = DockStyle.Fill;
            plotRiverDischarge.Location = new Point(3, 3);
            plotRiverDischarge.Margin = new Padding(4, 3, 4, 3);
            plotRiverDischarge.Name = "plotRiverDischarge";
            plotRiverDischarge.Size = new Size(1000, 727);
            plotRiverDischarge.TabIndex = 1;
            // 
            // tmr
            // 
            tmr.Interval = 2000;
            tmr.Tick += tmr_Tick;
            // 
            // butSave
            // 
            butSave.Location = new Point(12, 242);
            butSave.Name = "butSave";
            butSave.Size = new Size(129, 61);
            butSave.TabIndex = 3;
            butSave.Text = "Save";
            butSave.UseVisualStyleBackColor = true;
            butSave.Click += butSave_Click;
            // 
            // dlgSave
            // 
            dlgSave.FileName = "hoi";
            dlgSave.Filter = "csv|*.csv";
            // 
            // MainForm
            // 
            AutoScaleDimensions = new SizeF(7F, 15F);
            AutoScaleMode = AutoScaleMode.Font;
            ClientSize = new Size(1184, 761);
            Controls.Add(splitContainer1);
            Name = "MainForm";
            Text = "Data Optimization";
            splitContainer1.Panel1.ResumeLayout(false);
            splitContainer1.Panel2.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)splitContainer1).EndInit();
            splitContainer1.ResumeLayout(false);
            tabControl1.ResumeLayout(false);
            tabPage1.ResumeLayout(false);
            tabPage2.ResumeLayout(false);
            ResumeLayout(false);
        }

        #endregion

        private ScottPlot.FormsPlot plotRainfall;
        private SplitContainer splitContainer1;
        private Button butLoadData;
        private Button butRunOptimizer;
        private TabControl tabControl1;
        private TabPage tabPage1;
        private TabPage tabPage2;
        private ScottPlot.FormsPlot plotRiverDischarge;
        private Button butRunModel;
        private System.Windows.Forms.Timer tmr;
        private Button butSave;
        private SaveFileDialog dlgSave;
    }
}
