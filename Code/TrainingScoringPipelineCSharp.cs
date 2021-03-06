using System;
using System.Collections;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Remoting.Messaging;
using System.Runtime.Remoting.Metadata.W3cXsd2001;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using CommandLine;

namespace TrainingScoringPipeline
{
    class Program
    {
        class CmdlineOptions
        {
            [Option("trainingFile", Required = true, HelpText = "Training file")]
            public string trainingFile { get; set; }

            [Option("localTestFile", Required = true, HelpText = "Local test file")]
            public string localTestFile { get; set; }

            [Option("finalTestFile", Required = true, HelpText = "Final test file")]
            public string finalTestFile { get; set; }

            [Option("localTestTopRankerFile", Required = true, HelpText = "localTestTopRankerFile file")]
            public string localTestTopRankerFile { get; set; }

            [Option("localTestPredictionFile", Required = true, HelpText = "localTestPredictionFile file")]
            public string localTestPredictionFile { get; set; }

            [Option("finalTestPredictionFile", Required = true, HelpText = "finalTestPredictionFile file")]
            public string finalTestPredictionFile { get; set; }

            [Option("geneId", Required = true, HelpText = "geneId")]
            public string geneId { get; set; }

            [Option("paramSweepString", Required = true, HelpText = "paramSweepString")]
            public string paramSweepString { get; set; }

            [Option("paramSweepBatchSize", Required = true, HelpText = "paramSweepBatchSize")]
            public string paramSweepBatchSize { get; set; }

            [Option("paramSweepNumBatch", Required = true, HelpText = "paramSweepNumBatch")]
            public string paramSweepNumBatch { get; set; }

            [Option("bestModelZipFolder", Required = true, HelpText = "bestModelZipFolder")]
            public string bestModelZipFolder { get; set; }

            [Option("numGenesToTrain", Required = true, HelpText = "numGenesToTrain")]
            public string numGenesToTrain { get; set; }
        }



        //TODO
        //Arguments - --trainingFile ..\..\..\training.transposed.txt --localTestFile ..\..\..\Test\landmark.truth.merge.transpose.txt --finalTestFile ..\..\..\FinalTest\testData.transpose.txt --localTestTopRankerFile f1.txt --localTestPredictionFile f2.txt --finalTestPredictionFile f3.txt --geneId 971 --paramSweepString "{/p=lp{name=TREES min=500 max=1500 inc=100} /p=lp{name=LEAVES min=50 max=100 inc=20} /p=lp{name=MIL min=30 max=120 inc=20} /p=fp{name=SHRK min=0.25 max=4 logbase+} /p=fp{name=LR min=0.1 max=0.4 logbase+} /p=dp{name=TDROP v=0} s+ }" --paramSweepBatchSize 7 --paramSweepNumBatch 4 --bestModelZipFolder bestModelDir --numGenesToTrain 2

        // Test at 1, 10x, 100x stages while creating single modules.

        // Attach maml.exe, ResultProcessor.exe and dlls
        static void Main(string[] args)
        {
            // Get training file name, local test file name, final test file name, score debug output file, local gene score output, final test gene score output
            // Get the gene id and sweep parameters
            CmdlineOptions cmdLine = new CmdlineOptions();
            if (!CommandLine.Parser.Default.ParseArguments(args, cmdLine))
            {
                Console.WriteLine("Error in parsing command line options.");
                return;
            }

            int initGeneId = (int)new DataTable().Compute(cmdLine.geneId, null);
            cmdLine.geneId = initGeneId.ToString();

            //Total 12320 genes are present
            int highestGeneId = 12320;

            int genesToTrain = Int32.Parse(cmdLine.numGenesToTrain);

            for (int i = 0; i < genesToTrain; i++)
            {
                int curGeneId = initGeneId + i;
                Console.WriteLine("\n\n**Training Gene_" + curGeneId + "\n");

                if (curGeneId > highestGeneId)
                {
                    Console.WriteLine("\n\n**Request for GeneId > 12320. Request - Gene_" + curGeneId + "**\n\n");
                    continue;
                }

                string ModelTrainingOutputErrFolder = "SweepFolder_" + curGeneId + @"\";
                string ModelTrainingSummaryFile = ModelTrainingOutputErrFolder + "summary.txt";

                Directory.CreateDirectory("Folder_" + curGeneId);
                // Construct param sweep command
                string paramSweepCommand = "sweep sbs=" + cmdLine.paramSweepBatchSize + " snb=" +
                                           cmdLine.paramSweepNumBatch + " " +
                                           "sweeper=ldrand" + cmdLine.paramSweepString + " " +
                                           "runner=Local{ " +
                                           "pattern={TrainTest data=" + cmdLine.trainingFile + " " +
                                           "testfile=" + cmdLine.localTestFile + " " +
                                           "loader=TextLoader{sep=, col=Features:R4:0-969 col=Label:R4:" + (curGeneId - 1) +
                                           " header=+}" + " " +
                                           "tr=FastTreeRegression{lr=$LR$ shrk=$SHRK$ tdrop=$TDROP$ nl=$LEAVES$ iter=$TREES$ mil=$MIL$}" +
                                           " " +
                                           "sf=Folder_" + curGeneId + @"\summary_GENE" + curGeneId +
                                           ".out_TR$TREES$_LV$LEAVES$.txt " +
                                           "dout=Folder_" + curGeneId + @"\inst_GENE" + curGeneId +
                                           ".pred_TR$TREES$_LV$LEAVES$.txt " +
                                           "out=Folder_" + curGeneId + @"\model_GENE" + curGeneId +
                                           "_TR$TREES$_LV$LEAVES$.zip }" + " " +
                                           @"outfolder=" + ModelTrainingOutputErrFolder + " }";

                Console.WriteLine("\n\n**Constructecd sweep command\n" + paramSweepCommand + "\n\n");
                PrintDirInfo("Directory before paramsweep");

                DateTime start = DateTime.Now;

                // Execute sweep processor
                CreateMamlProcess(paramSweepCommand);

                DateTime end = DateTime.Now;
                Console.WriteLine("\n\n**Sweep Gene_" + curGeneId + " completed in: " + (end - start).Minutes + "min");

                PrintDirInfo("Directory after paramsweep and before resultsProcessor");
                Console.WriteLine("\n\n**Completed sweep command\n\n.");

                Console.WriteLine("**Start ResultProcessor**\n\n");

                // Execute ResultProcessor
                CreateResultProcessor(ModelTrainingOutputErrFolder, ModelTrainingSummaryFile);

                PrintDirInfo("Directory after resultprocessor");
                Console.WriteLine("\n\n**Completed ResultProcessor\n\n");

                DumpFullFileOnConsole(ModelTrainingSummaryFile);

                // Parse summary for top k ranker name and score.
                // Output top 5 best rankers. Add gene name column to data.
                // Ouput best ranker's score for local test data. Add gene name in header of single column
                string bestModelFilePath;
                ParseAndOutputTopKRankers(curGeneId, ModelTrainingSummaryFile, 5, cmdLine.localTestTopRankerFile,
                    cmdLine.localTestPredictionFile, out bestModelFilePath);

                string bestModelFileName = bestModelFilePath.Split('\\').Last(x => true);
                // Copy best model zip to output
                if (!Directory.Exists(cmdLine.bestModelZipFolder))
                {
                    Directory.CreateDirectory(cmdLine.bestModelZipFolder);
                }
                File.Copy(bestModelFilePath, cmdLine.bestModelZipFolder + @"\" + bestModelFileName);

                // Prepare scoring command
                string tempFinalPredictionFile = "Folder_" + curGeneId + @"\" + "finalprediction.temp.txt";
                string scoreRunCommand = "score " +
                                         "loader=TextLoader{sep=, col=Features:R4:0-969 header=+}" + " " +
                                         "data=" + cmdLine.finalTestFile + " " +
                                         "in=" + bestModelFilePath + " " +
                                         "dout=" + tempFinalPredictionFile + " ";

                Console.WriteLine("\n\n**Constructecd score command\n" + scoreRunCommand + "\n\n");
                PrintDirInfo("**Directory before score step**");

                // Execute score process for final test file
                CreateMamlProcess(scoreRunCommand);

                PrintDirInfo("**Directory after score step**");

                // Parse the gene score and output.
                ParseAndOutputFinalTestScore(curGeneId, tempFinalPredictionFile, cmdLine.finalTestPredictionFile);
            }

            // End execution
            Console.WriteLine("Finished Execution");
        }

        private static void ReturnWithGeneExceededMsg(int geneId, int highestGeneId, CmdlineOptions cmdLine)
        {
            string msg = "GeneId: " + geneId + " higher than " + highestGeneId + ". Not processed";
            Console.WriteLine(msg);
            WriteMsgToFile(msg, cmdLine.finalTestPredictionFile);
            WriteMsgToFile(msg, cmdLine.localTestPredictionFile);
            WriteMsgToFile(msg, cmdLine.localTestTopRankerFile);
            WriteMsgToFile(msg, cmdLine.bestModelZipFolder);
        }

        private static void WriteMsgToFile(string msg, string fileName)
        {
            using (StreamWriter sr = new StreamWriter(fileName))
            {
                sr.WriteLine(msg);
            }
        }

        private static void PrintDirInfo(string status)
        {
            Console.WriteLine("\n" + status);
            DirectoryInfo dir = new DirectoryInfo(Directory.GetCurrentDirectory());
            Console.WriteLine("Directory and file list - Hidding Dlls");
            foreach (DirectoryInfo d in dir.GetDirectories())
            {
                Console.WriteLine("{0, -30}\t directory", d.Name);
                foreach (var dc in d.GetFiles())
                {
                    Console.WriteLine("{0, -30}\t SubDirectoryFile", dc.Name);
                }
            }

            foreach (FileInfo f in dir.GetFiles())
            {
                if (!f.Name.EndsWith(".dll"))
                {
                    Console.WriteLine("{0, -30}\t File", f.Name);
                }

            }

            Console.WriteLine("\n");
        }

        private static void ParseAndOutputFinalTestScore(int geneId, string tempFinalPredictionFile, string finalPredictionFile)
        {
            /* #@ TextLoader{
                #@   header+
                #@   sep=tab
                #@   col=Score:R4:0
                #@ }
                Score
                4.390648
                4.290163
             * */

            using (StreamReader sr = new StreamReader(tempFinalPredictionFile))
            {
                StringReader predFileReader = null;
                if (File.Exists(finalPredictionFile))
                {
                    predFileReader = new StringReader(File.ReadAllText(finalPredictionFile));
                }

                using (StreamWriter sw = new StreamWriter(finalPredictionFile))
                {
                    string existingStr = "";
                    if (predFileReader != null)
                    {
                        existingStr = predFileReader.ReadLine();
                    }

                    string writeStr = existingStr + "," + "Gene_" + geneId;
                    sw.WriteLine(writeStr.Trim(new char[] { ',' }));

                    while (!sr.EndOfStream)
                    {
                        string line = sr.ReadLine();

                        if (line.Trim() == "")
                        {
                            continue;
                        }

                        double score;
                        if (double.TryParse(line, out score))
                        {
                            existingStr = "";
                            if (predFileReader != null)
                            {
                                existingStr = predFileReader.ReadLine();
                            }

                            writeStr = existingStr + "," + score.ToString();
                            sw.WriteLine(writeStr.Trim(new char[] { ',' }));
                        }
                    }
                }

            }
        }

        private static void DumpFullFileOnConsole(string filePath)
        {
            // Dump ModelTrainingSummaryFile
            using (StreamReader sr = new StreamReader(filePath))
            {
                Console.WriteLine("**Printing File**" + filePath);
                Console.WriteLine(sr.ReadToEnd());
            }
        }

        private static void ParseAndOutputTopKRankers(int geneId, string ModelTrainingSummaryFile, int topK, string topRankerOutputFile, string localTestBestPredictionFile, out string bestModelFilePath)
        {
            using (StreamReader sr = new StreamReader(ModelTrainingSummaryFile))
            {
                // Skip 2 lines
                sr.ReadLine();
                string summaryheader = sr.ReadLine();

                int l2ErrorIdx = summaryheader.Split('\t').ToList().FindIndex(x => x.Contains("L2(avg)"));
                int commandIdx = summaryheader.Split('\t').ToList().FindIndex(x => x.Contains("Command"));
                Console.WriteLine("In summary file L2(Avg) found at idx {0}, command line found at idx {1}", l2ErrorIdx, commandIdx);

                List<Tuple<double, string>> l2ErrorList = new List<Tuple<double, string>>();
                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine();

                    if (line.Trim() == "")
                    {
                        continue;
                    }

                    string[] lineparts = line.Split('\t');
                    l2ErrorList.Add(new Tuple<double, string>(Double.Parse(lineparts[l2ErrorIdx]), line));
                }

                // lowest error at start
                l2ErrorList.Sort((x1, x2) => x1.Item1.CompareTo(x2.Item1));

                // Print best model
                Console.WriteLine("\n\n**Top model Gene_{0}\n{1}", geneId, l2ErrorList[0].Item2);

                // lowest error command is here. Get best model file name
                string command = l2ErrorList[0].Item2.Split('\t')[commandIdx];
                string[] splitcommand = command.Split('=');
                bestModelFilePath = splitcommand[splitcommand.Length - 1].Trim();

                // write topk models
                using (StreamWriter sw1 = new StreamWriter(topRankerOutputFile, append:true))
                {
                    sw1.WriteLine("GeneId\t" + summaryheader);
                    for (int i = 0; i < topK && i < l2ErrorList.Count; i++)
                    {
                        sw1.WriteLine("Gene_" + geneId + "\t" + l2ErrorList[i].Item2);
                    }
                }

                // write the local test prediction values
                string instOutputFile = command.Split(' ').First(x => x.StartsWith("dout="));
                instOutputFile = instOutputFile.Split('=')[1];

                using (StreamReader sr1 = new StreamReader(instOutputFile))
                {
                    // skip header
                    sr1.ReadLine();

                    StringReader existingPredFileReader = null;
                    if (File.Exists(localTestBestPredictionFile))
                    {
                        existingPredFileReader = new StringReader(File.ReadAllText(localTestBestPredictionFile));
                    }

                    using (StreamWriter sw2 = new StreamWriter(localTestBestPredictionFile))
                    {
                        String existingString = "";
                        if (existingPredFileReader != null)
                        {
                            existingString = existingPredFileReader.ReadLine();
                        }

                        string writeStr = existingString + "," + "Gene_" + geneId;
                        sw2.WriteLine(writeStr.Trim(new char[]{','}));

                        while (!sr1.EndOfStream)
                        {
                            string line = sr1.ReadLine();
                            if (line.Trim() == "")
                            {
                                continue;
                            }

                            existingString = "";
                            if (existingPredFileReader != null)
                            {
                                existingString = existingPredFileReader.ReadLine();
                            }

                            writeStr = existingString + "," + line.Split('\t')[2];
                            sw2.WriteLine(writeStr.Trim(new char[]{','}));
                        }
                    }
                }
            }
        }

        private static void CreateResultProcessor(string ModelTrainingOutputErrFolder, string ModelTrainingSummaryFile)
        {
            Process p = new Process();

            ProcessStartInfo SI = new ProcessStartInfo();
            SI.CreateNoWindow = false;
            SI.FileName = "ResultProcessor.exe";
            SI.UseShellExecute = false;
            SI.Arguments = ModelTrainingOutputErrFolder + "*.out.txt outputFile=" + ModelTrainingSummaryFile;

            p.StartInfo = SI;
            try
            {
                p.Start();
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
            }

            p.WaitForExit();

        }

        private static void CreateMamlProcess(string command)
        {
            Process p = new Process();

            ProcessStartInfo SI = new ProcessStartInfo();
            SI.CreateNoWindow = false;
            SI.FileName = "maml.exe";
            SI.UseShellExecute = false;
            SI.Arguments = command;

            p.StartInfo = SI;
            try
            {
                p.Start();
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
            }

            p.WaitForExit();
        }

        // Create a column joiner module.
    }
}
