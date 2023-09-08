#include <iostream>
#include <string>
#include <map>
#include <fstream> 
#include <nlohmann/json.hpp>
#include <cmath> 
#include <set>

void establishEvents(std::ifstream& reads_bed_file, std::map<std::string, std::map<int, double>>& intervals, const std::string& analysis_type, const int& userStart, const int& userEnd, double& correctionFactor, std::vector<std::string>& geneNames, std::set<std::string>& uniqueGeneNameChecker);
void computeCounts(std::map<std::string, std::map<int, double>>& intervals, std::map<std::string, std::map<int, double>>& intervalsCounts);
void fillInMissingPositions(std::map<std::string, std::map<int, double>>& intervalsCounts, std::map<std::string, std::map<int, double>>& intervalsCountsFinal, const int& userEnd);
void WeightedAverageSTDEV(std::map<std::string, std::map<int, double>>& intervalsCountsFinal, std::map<std::string, std::array<double, 3>>& geneWightedAvgStdev);
void outputTheFile(std::ofstream& outputFile, std::map<std::string, std::array<double, 3>>& geneWightedAvgStdev, const std::string& analysis_type, std::vector<std::string>& geneNames, const std::string& outputName); 

int main(int argc, char* argv[]) {
    // make sure that 6 arguments are passed
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <json data> <type of analysis(5', 3' or centers)> <output name> <start of regions coordinate (-500)> <end of regions coordinate (500)>" << std::endl << std::endl;
        exit(EXIT_FAILURE);
        return 1;
    }
    // check for valid type_of_overlap option and establish type of analysis; declaring bedFilePath to be opened, jsonData to be parsed 
    std::string receivedJsonString = argv[1];
    std::string analysis_type = argv[2];
    std::string outputName = argv[3];
    int userStart = std::stoi(argv[4]);
    int userEnd = std::stoi(argv[5]);

    nlohmann::json data = nlohmann::json::parse(receivedJsonString);

    // create vector of Interval objects
    std::map<std::string, std::map<int, double>> intervals;
    // vector to store gene names in order
    std::vector<std::string> geneNames;
    // set to store unique gene names
    std::set<std::string> uniqueGeneNameChecker;
    try {
        for (auto& element : data.items()) {
            std::string bedFilePath = element.key();
            
            // Check if the value is a number
            if (element.value().is_number()) {
                double correctionFactor = element.value().get<double>();
                std::ifstream reads_bed_file(bedFilePath);
                if (!reads_bed_file.is_open()) {
                    std::cerr << "Unable to open file: " << bedFilePath << std::endl;
                    continue;
                }
                establishEvents(reads_bed_file, intervals, analysis_type, userStart, userEnd, correctionFactor, geneNames, uniqueGeneNameChecker);

                reads_bed_file.close();
            } else {
                std::cerr << "Unexpected value for key '" << bedFilePath << "'. Expected a number." << std::endl;
            }
        }
        } catch (const nlohmann::json::exception& e) {
            std::cerr << "JSON error: " << e.what() << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Standard exception: " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknown error occurred!" << std::endl;
        }

    
    std::map<std::string, std::map<int, double>> intervalsCounts;
    computeCounts(intervals, intervalsCounts); 
    // done with intervals, so clear it
    intervals.clear();
    // make final map with all positions filled in
    std::map<std::string, std::map<int, double>> intervalsCountsFinal;  
    fillInMissingPositions(intervalsCounts, intervalsCountsFinal, userEnd);
    // done with intervalsCounts, so clear it
    intervalsCounts.clear();
    // calculate weighted average
    std::map<std::string, std::array<double, 3>> geneWightedAvgStdev;
    WeightedAverageSTDEV(intervalsCountsFinal, geneWightedAvgStdev);
    // output file
    std::ofstream outputFile;
    outputTheFile(outputFile, geneWightedAvgStdev, analysis_type, geneNames, outputName);

    return 0;
}


void establishEvents(std::ifstream& reads_bed_file, std::map<std::string, std::map<int, double>>& intervals, const std::string& analysis_type, const int& userStart, const int& userEnd, double& correctionFactor, std::vector<std::string>& geneNames, std::set<std::string>& uniqueGeneNameChecker){
    // read each line of the reads bed file
    // declare variables to store bed file information from bedtools output
    std::string regionChr; // col1
    int regionStart; // col2
    int regionEnd; // col3
    std::string regionGene; // col4
    std::string additionalInfo1; // col5
    std::string regionStrand; // col6      
    std::string readsChr; // col7
    int readsStart; // col8
    int readsEnd; // col9
    std::string readsID; // col 10 
    int readsQual; // col 11
    std::string readsStrand; // col12 

    while (reads_bed_file >> regionChr >> regionStart >> regionEnd >> regionGene >> additionalInfo1 >> regionStrand >> readsChr >> readsStart >> readsEnd >> readsID >> readsQual >> readsStrand) {
        // add genes to vector
        if (uniqueGeneNameChecker.insert(regionGene).second) {
            // .second of the pair returned by .insert() is true if the insertion took place, i.e., the item was unique
            geneNames.push_back(regionGene);
        }
        int readBasePositionStart;
        int readBasePositionEnd;
        int positionTSS = (regionStart + regionEnd) / 2;
        float readCenter;
        // calculate read position relative to TSS in the middle of the region
        if (analysis_type == "5"){
            // If positive strand, start is the 5' end, end is the 5' end + 1
             // If negative strand, start is the 3' end - 1, end is the 3' end 

            if (regionStrand == "+"){
                readBasePositionStart = readsStart - positionTSS;
                readBasePositionEnd = readBasePositionStart;
            }

            else if (regionStrand == "-"){
                readBasePositionStart = positionTSS - readsEnd;
                readBasePositionEnd = readBasePositionStart;
            }

        }
        else if (analysis_type == "3"){
            // If positive strand, start is the 3' end - 1, end is the 3' end
            // If negative strand, start is the 5' end, end is the 5' end + 1

            if (regionStrand == "+"){
                readBasePositionStart = readsEnd - positionTSS;
                readBasePositionEnd = readBasePositionStart;
            }

            else if (regionStrand == "-"){
                readBasePositionStart = positionTSS - readsStart;
                readBasePositionEnd = readBasePositionStart;
            }
        }

        else if (analysis_type == "centers"){
            if (regionStrand == "+"){
                readBasePositionStart = ((readsStart + readsEnd) / 2) - positionTSS;  
                readBasePositionEnd = readBasePositionStart;
            }

            else if (regionStrand == "-"){
                readBasePositionStart = positionTSS - (readsStart + readsEnd) / 2; 
                readBasePositionEnd = readBasePositionStart;
            }
        }

        else if (analysis_type == "full"){
            if (regionStrand == "+"){
                readBasePositionStart = readsStart - positionTSS;
                readBasePositionEnd = readsEnd - positionTSS;
            }

            else if (regionStrand == "-"){
                readBasePositionStart = positionTSS - readsEnd;
                readBasePositionEnd = positionTSS - readsStart;
            }
        }

        // if the postion is 0 or positive, add 1 to the position to reflect starting at +1
        if (readBasePositionStart >= 0){
            readBasePositionStart += 1;
        }
        if (readBasePositionEnd >= 0){
            readBasePositionEnd += 1;
        }

        // now filter out reads that are outside of the user defined start and end
        // 4 types of reads to count:
        int newStart = -99999;
        int newEnd = -99999;
        // 1. both left and right position of read are outside of user defined start and end but read overlaps with user defined start and end
        if (readBasePositionStart < userStart && readBasePositionEnd > userEnd){
            // now add the read to the vector of reads for that interval
            newStart = userStart;
            newEnd = userEnd;
        }
        // 2. left position of read is less than start but right position of read is inside user defined start and end
        else if (readBasePositionStart < userStart && readBasePositionEnd >= userStart && readBasePositionEnd <= userEnd){
            newStart = userStart;
            newEnd = readBasePositionEnd;
        }
        // 3. left position of read is inside user defined start and end but right position of read is greater than end
        else if (readBasePositionStart >= userStart && readBasePositionStart <= userEnd && readBasePositionEnd > userEnd){
            newStart = readBasePositionStart;
            newEnd = userEnd;
        }
        // 4. both left and right position of read are inside user defined start and end
        else if (readBasePositionStart >= userStart && readBasePositionStart <= userEnd && readBasePositionEnd >= userStart && readBasePositionEnd <= userEnd){
            newStart = readBasePositionStart;
            newEnd = readBasePositionEnd;
        }

        if (newStart != -99999 && newEnd != -99999){
            // addto intervals
            if (intervals[regionGene].find(newStart) == intervals[regionGene].end()) {
                    intervals[regionGene][newStart] = 1 * correctionFactor; // add 1 to the count and normalize by the correction factor

            } else {
                intervals[regionGene][newStart] += 1 * correctionFactor; 
            }

            int newEndNotZero; // for cases where newStart is -1, newEnd is 0, we want to add 1 to newEnd
            if ((newEnd+1) == 0){
                newEndNotZero = 1;
            }
            else{
                newEndNotZero = newEnd+1;
            }
            if (intervals[regionGene].find(newEndNotZero) == intervals[regionGene].end()) {

                intervals[regionGene][newEndNotZero] = -1 * correctionFactor;
            } else {
                intervals[regionGene][newEndNotZero] -= 1 * correctionFactor;
            }
        }
    }
}

void computeCounts(std::map<std::string, std::map<int, double>>& intervals, std::map<std::string, std::map<int, double>>& intervalsCounts) {
    // iterate through each key-value pair in the map
    for (const auto& [gene, intervalPairs]: intervals){
        // initialize current count
        double currentCount = 0;
        // iterate through each pair in the vector
        for (const auto& [position, change]: intervalPairs){
            // update the current count
            currentCount += change;
            // update the result in intervalsCounts
            intervalsCounts[gene][position] = currentCount;
        }
    }
}

void fillInMissingPositions(std::map<std::string, std::map<int, double>>& intervalsCounts, std::map<std::string, std::map<int, double>>& intervalsCountsFinal, const int& userEnd){
    // iterate through each key-value pair in the map
    for (const auto& [gene, intervalPairs]: intervalsCounts) {
        // initialize vector with seen positions
        std::vector<int> seenPositions;
        // first loop: populate seenPositions
        for (const auto& [position, count]: intervalPairs) {
            seenPositions.push_back(position);
        }
        // second loop: iterate through seenPositions to fill missing positions
        for (size_t eachPosition = 0; eachPosition < seenPositions.size() - 1; eachPosition++) {
            for (int absentPositions = seenPositions[eachPosition] + 1; 
                 absentPositions < seenPositions[eachPosition + 1]; absentPositions++) {
                if (absentPositions != 0) {  // skip if absentPositions is 0
                    intervalsCountsFinal[gene][absentPositions] = intervalsCounts[gene][seenPositions[eachPosition]];
                }
            }
        }

        // copy existing positions to the final map
        for (const auto& [position, count] : intervalPairs) {
            if (position <= userEnd){ // skip last position if is greater than userEnd
                intervalsCountsFinal[gene][position] = count;
            }
        }
    }
}

void WeightedAverageSTDEV(std::map<std::string, std::map<int, double>>& intervalsCountsFinal, std::map<std::string, std::array<double, 3>>& geneWightedAvgStdev){
    for (const auto& [gene, values] : intervalsCountsFinal) {
        double sum = 0;
        double totalWeights = 0;

        for (const auto& [base, count] : values) {
            sum += base * count;
            totalWeights += count;
        }

        double weightedAverage = (totalWeights != 0) ? sum / totalWeights : 0;

        double varianceSum = 0;
        for (const auto& [base, count] : values) {
            varianceSum += count * (base - weightedAverage) * (base - weightedAverage);
        }

        double variance = (totalWeights != 0) ? varianceSum / totalWeights : 0;
        double stdev = std::sqrt(variance);
        
        geneWightedAvgStdev[gene] = {static_cast<double>(totalWeights), weightedAverage, stdev};
    }
}

void outputTheFile(std::ofstream& outputFile, std::map<std::string, std::array<double, 3>>& geneWightedAvgStdev, const std::string& analysis_type, std::vector<std::string>& geneNames, const std::string& outputName){
    outputFile.open(outputName + "wighted_average_stdev.bed");
    outputFile << "Gene\tCounts\tWeighted-average\tSTDEV" << std::endl;
    for (const auto& gene : geneNames){
        // if gene not in the map, add it with 0 values
        if (geneWightedAvgStdev.find(gene) == geneWightedAvgStdev.end()) {
            geneWightedAvgStdev[gene] = {0, 0, 0};
        }
        outputFile << gene << "\t" << geneWightedAvgStdev[gene][0] << "\t" << geneWightedAvgStdev[gene][1] << "\t" << geneWightedAvgStdev[gene][2] << std::endl;
    }
    outputFile.close();
}
