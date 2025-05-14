#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "ns3/lorawan-module.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <queue>
#include <ctime>
#include <utility>
#include <string>
#include <iomanip>
#include <set>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <limits>


using namespace ns3;
using namespace ns3::lorawan;
using namespace std;

int numNodes = 10;
int areaSize = 10000;
double pursuit = 200.0;
double optimalRange = 1000.0;
double shrinkFactor = 300.0;
int logfull;
int logSteps;
int stage;
int minDistConst = areaSize * 2.0;
vector<vector<Time>> packetTrackers;
vector<vector<vector<double>>> dataRates;
vector<Ptr<NetDevice>> devicePtrs;
vector<Ptr<MobilityModel>> models;
vector<set<int>> links;
LoraPhyHelper phyHelper;
LorawanMacHelper macHelper;
LoraHelper lora;
NodeContainer nodes;


double Distance(const Vector& a, const Vector& b)
{
    return std::sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

vector<pair<int, int>> Subgraphs(const vector<set<int>>& graph, const vector<Ptr<MobilityModel>>& models)
{
    int n = graph.size();
    vector<bool> visited(n, false);
    vector<set<int>> subgraphs;

    // Find connected components via DFS
    function<void(int, set<int>&)> dfs = [&](int node, set<int>& component) {
        visited[node] = true;
        component.insert(node);
        for (int neighbor : graph[node]) {
            if (!visited[neighbor]) {
                dfs(neighbor, component);
            }
        }
    };

    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            set<int> component;
            dfs(i, component);
            subgraphs.push_back(component);
        }
    }

    if (logSteps != 0) std::cout << "Subgraphs count: " << subgraphs.size() << "\n";

    // Find closest pairs between all unique pairs of subgraphs
    vector<pair<int, int>> closestPairs;
    for (size_t i = 0; i < subgraphs.size(); ++i) {
        for (size_t j = i + 1; j < subgraphs.size(); ++j) {
            double minDist = minDistConst;
            pair<int, int> closest;

            for (int a : subgraphs[i]) {
                for (int b : subgraphs[j]) {
                    int dist = Distance(models[a]->GetPosition(), models[b]->GetPosition());
                    if (dist < minDist) {
                        minDist = dist;
                        closest = {a, b};
                    }
                }
            }

            closestPairs.push_back(closest);
        }
    }

    return closestPairs;
}




void DrawAsciiMap(vector<Ptr<MobilityModel>> positions, double areaSize) {
    const int mapWidth = 50;
    const int mapHeight = 50;
    const double cellWidth = areaSize / mapWidth;
    const double cellHeight = areaSize / mapHeight;

    vector<vector<char>> grid(mapHeight, vector<char>(mapWidth, '.'));

    for (size_t i = 0; i < positions.size() - 1; ++i) {
        Vector pos = positions[i]->GetPosition();
        int col = min(static_cast<int>(pos.x / cellWidth), mapWidth - 1);
        int row = min(static_cast<int>(pos.y / cellHeight), mapHeight - 1);
        grid[mapHeight - 1 - row][col] = '0' + i;
    }
    Vector pos = positions[positions.size() - 1]->GetPosition();
    int col = min(static_cast<int>(pos.x / cellWidth), mapWidth - 1);
    int row = min(static_cast<int>(pos.y / cellHeight), mapHeight - 1);
    grid[mapHeight - 1 - row][col] = 'X';

    std::cout << "ASCII Map (" << areaSize << " x " << areaSize << " meters):\n";
    for (const auto& row : grid) {
        for (char cell : row) {
            std::cout << cell << ' ';
        }
        std::cout << "\n";
    }
}


double Dot(const Vector& v1, const Vector& v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

int RestorationAlgorithm(vector<Ptr<MobilityModel>> mobilityNodes, double areaSize, int samplingRes = 100)
{
    set<int> primaryAffected;
    for (int i = 0; i < numNodes; i++)
    {
        int connections = 0;
        for (int j = 0; j < numNodes; j++)
        {
            if (dataRates[stage][j][i] != 0 && dataRates[stage][i][j] != 0) connections++;
        }
        if (connections < 2) primaryAffected.insert(i);
    }

    int result = 0;
    vector<set<int>> newLinks;
    newLinks.resize(numNodes);
    for (int i = 0; i < dataRates[stage].size(); i++) {
	for (int j = i + 1; j < dataRates[stage][i].size(); j++)
	{
            if (dataRates[stage][i][j] != 0.0 && dataRates[stage][j][i] != 0.0)
            {
                newLinks[i].insert(j);
                newLinks[j].insert(i);
            }
        }
    }

    if (primaryAffected.size() > 0)
    {
        result = primaryAffected.size();
        vector<vector<Vector>> regions(numNodes + 1); // Voronoi regions
        Vector ifPos = mobilityNodes[numNodes]->GetPosition();

        // Sample the space and assign points to closest node
        for (int i = 0; i < samplingRes; ++i)
	{
            for (int j = 0; j < samplingRes; ++j)
	    {
                Vector sample = { (1.0 * i * areaSize) / samplingRes, (1.0 * j * areaSize) / samplingRes, 0 };
                int closest = numNodes; // closest index set to interference source
                Vector closestPos = ifPos; // closest point set to interference source
                double minDist = Distance(sample, closestPos);

                for (int n = 0; n < numNodes; ++n)
		{
                    Vector pos = mobilityNodes[n]->GetPosition();
                    double d = Distance(sample, pos);
                    if (d < minDist)
		    {
                        minDist = d;
                        closest = n;
                    }
                }
                regions[closest].push_back(sample);
            }
        }

        // Move each node to the centroid of its region
        for (int n = 0; n < numNodes; ++n)
	{
            if (regions[n].empty()) continue;

            double sumX = 0, sumY = 0;
            for (const Vector& p : regions[n])
	    {
                sumX += p.x;
                sumY += p.y;
            }
            Vector pointPos = mobilityNodes[n]->GetPosition();
            Vector lloyds = {sumX / regions[n].size(), sumY / regions[n].size(), 0};
            double dot = Dot(lloyds - pointPos, ifPos - pointPos) / Distance(lloyds, pointPos) / Distance(ifPos, pointPos);

            if (primaryAffected.count(n))
            {
                Vector correctionA = {0, 0, 0};
                if (dot >= -0.5) correctionA = (pointPos - ifPos) * ( pursuit / Distance(pointPos, ifPos)); // check
                else correctionA = (lloyds - pointPos) * ( pursuit / Distance(pointPos, lloyds));

                Vector mean = {0, 0, 0};
                double minDist = minDistConst;
                double minDistConnected = minDistConst;
                int closest = -1;
                int closestConnected = -1;
                for(int peer : links[n])
                {
                    Vector peerPos = mobilityNodes[peer]->GetPosition();
                    mean = mean + peerPos;
                    if (Distance(peerPos, pointPos) < minDist && !newLinks[n].count(peer) && !primaryAffected.count(peer))
                    {
                        closest = peer;
                        minDist = Distance(peerPos, pointPos);
                    }
                    if (Distance(peerPos, pointPos) < minDistConnected && newLinks[n].count(peer))
                    {
                        closestConnected = peer;
                        minDistConnected = Distance(peerPos, pointPos);
                    }
                }
                if (logSteps != 0) std::cout << n <<"'s closest disconnected peer: " << closest << "\n";
                if (logSteps != 0) std::cout << n <<"'s closest connected peer: " << closestConnected << "\n";
                if (closestConnected != -1 && minDistConnected < optimalRange)
                {
                    Vector closestConnectedPos = mobilityNodes[closestConnected]->GetPosition();
                    mobilityNodes[closestConnected]->SetPosition(closestConnectedPos - (pointPos - closestConnectedPos) * ((optimalRange - minDistConnected) / 2.0 / Distance(closestConnectedPos, pointPos))); //check
                }
                if (closest != -1 && minDist > 2 * shrinkFactor)
                {
                    Vector closestPos = mobilityNodes[closest]->GetPosition();
                    mobilityNodes[closest]->SetPosition(closestPos + (pointPos - closestPos) * (shrinkFactor / Distance(closestPos, pointPos))); // check
                }
                mean = mean * (1.0 / links[n].size());
                Vector correctionB = (mean - pointPos) * (pursuit / Distance(mean, pointPos)); // check
                double corrDot = Dot(correctionA, correctionB);

                if (corrDot < 0) correctionA = correctionA - ((corrDot / correctionB.GetLength() / correctionB.GetLength()) * correctionB);

                Vector final = pointPos + correctionA + correctionB;
                final.x = std::max(0.0, final.x);
                final.x = std::min(areaSize, final.x);
                final.y = std::max(0.0, final.y);
                final.y = std::min(areaSize, final.y);

                mobilityNodes[n]->SetPosition(final);
            }
        }
    }
    else
    {
        vector<pair<int,int>> pairs = Subgraphs(newLinks, mobilityNodes);
        result = pairs.size();
        for (pair<int,int> couple : pairs)
        {
            Vector first = mobilityNodes[couple.first]->GetPosition();
            Vector second = mobilityNodes[couple.second]->GetPosition();
            double dist = Distance(first, second);
            mobilityNodes[couple.first]->SetPosition(first + (second - first) * (dist / 4.0 / Distance(first, second)));
            mobilityNodes[couple.second]->SetPosition(second + (first - second) * (dist / 4.0 / Distance(first, second)));
        }
    }
    return result;
}


// Packs two non-negative integers into a single unique integer
unsigned long long cantorPair(unsigned long long a, unsigned long long b)
{
    return ((a + b) * (a + b + 1)) / 2 + b;
}

// Unpacks the single integer back into the original two integers
pair<unsigned long long, unsigned long long> cantorUnpair(unsigned long long z)
{
    unsigned long long w = static_cast<unsigned long long>((sqrt(8.0 * z + 1) - 1) / 2);
    unsigned long long t = (w * (w + 1)) / 2;
    unsigned long long b = z - t;
    unsigned long long a = w - b;
    return {a, b};
}

void OnPhyTxStart(Ptr<const Packet> packet, uint32_t sf)
{
    FlowIdTag tag;
    if (packet->PeekPacketTag(tag))
    {
        int flow = tag.GetFlowId();
        int src = cantorUnpair(flow).first;
        int dst = cantorUnpair(flow).second;
        packetTrackers[src][dst] = Simulator::Now();
    }
}

void OnPhyPacketReceived(Ptr<NetDevice> device, Ptr<const Packet> packet, uint32_t sf)
{
    FlowIdTag tag;
    int deviceIdx;
    for (size_t i = 0; i < devicePtrs.size(); ++i)
    {
        if (devicePtrs[i] == device)
        {
            deviceIdx = i;
            break;
        }
    }
    if (packet->PeekPacketTag(tag))
    {
        int flow = tag.GetFlowId();
        int expectedSrc = cantorUnpair(flow).first;
        int expectedDst = cantorUnpair(flow).second;
        if (expectedDst == deviceIdx)
        {
            Time start = packetTrackers[expectedSrc][expectedDst];
            Time end = Simulator::Now();
            double duration = (end - start).GetSeconds();
            double bps = (packet->GetSize() * 8.0) / duration;
            dataRates[stage][expectedSrc][expectedDst] = bps;
            if (logfull != 0) std::cout << "Node " << expectedDst << " received packet from " << expectedSrc << ", time : " << duration << ", bps: " << bps << "\n";
        }
        else
        {
            if (logfull != 0) std::cout << "Node " << deviceIdx << " received unwanted packet from " << expectedSrc << "\n";
        }
    }
}

double Distance(Ptr<MobilityModel> a, Ptr<MobilityModel> b)
{
    return a->GetDistanceFrom(b);
}

void ApplyInterference(int startDelay, int lengthTenthSeconds, double power, int sf)
{
    if (logfull != 0) std::cout << "Interference will start at " << startDelay << "\n";
    for (int i = 0; i < lengthTenthSeconds; i++)
    {
        Simulator::Schedule(Seconds(startDelay + i * 0.1), [=]() {
            Ptr<Packet> ifPacket = Create<Packet>(100);
            Ptr<LoraNetDevice> loraInterferenceDev = DynamicCast<LoraNetDevice>(devicePtrs[numNodes]);

            double frequency = 868100000;
            LoraTxParameters txIfParams;
            txIfParams.sf = sf;
            double txIfPower = power;
            loraInterferenceDev->GetPhy()->Send(ifPacket, txIfParams, frequency, txIfPower);
        });
    }
    if (logfull != 0) std::cout << "Interference will stop at " << startDelay + lengthTenthSeconds * 0.1 << "\n";
}

int NetworkCheck(int startDelay)
{
    int oneSimDuration = 10;
    int schedule = -1;
    for (int src = 0; src < numNodes; src++) {
        for (int dst = 0; dst < numNodes; dst++) {
            if (src == dst) continue;

            for (int sf = 12; sf >= 7; sf--)
            {
                schedule++;
                int scheduleOffset = startDelay + schedule * oneSimDuration;
                if (logfull != 0) std::cout << "Scheduled at " << scheduleOffset << " sec, Source: " << src << ", Destination: " << dst << ", SF: " << sf << "\n";
                Simulator::Schedule(Seconds(scheduleOffset), [=]() {
                    Ptr<Packet> packet = Create<Packet>(100);
                    Ptr<LoraNetDevice> loraDev = DynamicCast<LoraNetDevice>(devicePtrs[src]);
                    Ptr<Packet> ifPacket = Create<Packet>(100);
                    Ptr<LoraNetDevice> loraInterferenceDev = DynamicCast<LoraNetDevice>(devicePtrs[numNodes]);

                    FlowIdTag flowTag;
                    flowTag.SetFlowId(cantorPair(src,dst));
                    packet->AddPacketTag(flowTag);

                    LoraTxParameters txParams;
                    double dist = Distance(models[src], models[dst]);
                    txParams.sf = sf;
                    txParams.bandwidthHz = 125000;
                    double frequency = 868100000;
                    double txPower = 14.0;
                    if (logfull != 0) std::cout << "Sending packet: " << src << " -> " << dst << " with SF = " << sf << ", and distance: " << dist << "\n";
                    loraDev->GetPhy()->Send(packet, txParams, frequency, txPower);
                });
            }
        }
    }
    return (schedule + 1) * 10;
}

void PrintData()
{
    std::cout << "\nNode Positions:\n";
    for (int i = 0; i < numNodes; ++i)
    {
        Vector pos = models[i]->GetPosition();
        std::cout << "Node " << i << ": (" << pos.x << ", " << pos.y << ")\n";
    }
    Vector pos = models[numNodes]->GetPosition();
    std::cout << "Interference node: ( " << pos.x << ", " << pos.y << ")\n";

    std::cout << "\nMesh Connections:\n";
    for (int i = 0; i < dataRates[stage].size(); i ++)
    {
        std::cout << "Node " << i << " connects to: ";
        for (int j = 0; j < dataRates[stage][i].size(); j++)
	{
            if (dataRates[stage][i][j] != 0.0 && dataRates[stage][j][i] != 0.0) std::cout << j << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nData Rate Matrix (bps):\n    ";
    int cellWidth = 10;

    for (int j = 0; j < numNodes; ++j)
        std::cout << std::setw(cellWidth) << j;
    std::cout << "\n";

    for (int i = 0; i < numNodes; ++i)
    {
        std::cout << std::setw(3) << i << ":";
        for (int j = 0; j < numNodes; ++j)
	{
            std::cout << std::setw(cellWidth) << dataRates[stage][i][j];
        }
        std::cout << "\n";
    }
    DrawAsciiMap(models, areaSize);
    std::cout << "\n";
}

void Setup()
{
    NetDeviceContainer devices = lora.Install(phyHelper, macHelper, nodes);

    for (uint32_t i = 0; i < devices.GetN(); ++i) {
        devicePtrs.push_back(devices.Get(i));
    }

    for (int i = 0; i < numNodes; ++i) {
        Ptr<LoraNetDevice> dev = DynamicCast<LoraNetDevice>(devicePtrs[i]);
        dev->GetPhy()->TraceConnectWithoutContext("StartSending", MakeCallback(&OnPhyTxStart));
        dev->GetPhy()->TraceConnectWithoutContext("ReceivedPacket", MakeBoundCallback(&OnPhyPacketReceived, dev));
    }
}

int main(int argc, char *argv[]) {
    CommandLine cmd;
    int rand = 0;
    double power = 0.0;
    int sf = 0;
    logfull = 0;
    logSteps = 0;
    stage = 0;
    int maxIter = 100;
    numNodes = 10;
    areaSize = 10000;

    cmd.AddValue("power", "Interference power", power);
    cmd.AddValue("sf", "Interference spreading factor", sf);
    cmd.AddValue("rand", "Random", rand);
    cmd.AddValue("logfull", "Full logs", logfull);
    cmd.AddValue("logsteps", "Steps logs", logSteps);
    cmd.Parse(argc, argv);

    if (rand == 0) SeedManager::SetSeed(time(NULL));
    else SeedManager::SetSeed(rand);
    std::cout << "Rand seed: " << time(NULL) << "\n";
    SeedManager::SetRun(1);

    packetTrackers.resize(numNodes, std::vector<Time>(numNodes, Simulator::Now()));

    nodes.Create(numNodes + 1);
    dataRates.resize(maxIter + 1, std::vector<std::vector<double>>(numNodes, std::vector<double>(numNodes, 0.0)));
    links.resize(numNodes);

    std::string xRange = "ns3::UniformRandomVariable[Min=0.0|Max=" + std::to_string(areaSize) + "]";
    std::string yRange = "ns3::UniformRandomVariable[Min=0.0|Max=" + std::to_string(areaSize) + "]";

    MobilityHelper mobility;
    mobility.SetPositionAllocator("ns3::RandomRectanglePositionAllocator", "X", StringValue(xRange), "Y", StringValue(yRange));
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    mobility.Install(nodes);

    models.resize(numNodes + 1);
    for (int i = 0; i < numNodes + 1; ++i) {
        models[i] = nodes.Get(i)->GetObject<MobilityModel>();
    }
    Vector ifPos = models[numNodes]->GetPosition();
    Vector adjustment = {areaSize / 4, areaSize / 4, 0};
    models[numNodes]->SetPosition(ifPos * 0.5 + adjustment);

    Ptr<LogDistancePropagationLossModel> loss = CreateObject<LogDistancePropagationLossModel>();
    Ptr<ConstantSpeedPropagationDelayModel> delay = CreateObject<ConstantSpeedPropagationDelayModel>();
    Ptr<LoraChannel> channel = CreateObject<LoraChannel>(loss, delay);

    phyHelper.SetChannel(channel);
    phyHelper.SetDeviceType(LoraPhyHelper::GW);
    macHelper.SetDeviceType(LorawanMacHelper::GW);

    Setup();

    int baseDelay = 1;
    int duration = NetworkCheck(baseDelay);
    Simulator::Run();
    PrintData();

    for (int i = 0; i < dataRates[0].size(); i++) {
        for (int j = i + 1; j < dataRates[0][i].size(); j++) {
            if (dataRates[0][i][j] != 0.0 && dataRates[0][j][i] != 0.0)
            {
                links[i].insert(j);
                links[j].insert(i);
            }
        }
    }
    for (int i = 0; i < links.size(); i++)
    {
        if (links[i].size() < 2) return 0;
    }
    if (Subgraphs(links, models).size() > 0) return 0;

    int affectedCount = 1;
    for (int i = 0; i < maxIter && affectedCount > 0; i++)
    {
        stage++;
        duration = NetworkCheck(baseDelay);
        ApplyInterference(0, (duration + 2) * 10, power, sf);
        Simulator::Run();
        if (i == 0) PrintData();
        if (logSteps != 0) PrintData();
        affectedCount = RestorationAlgorithm(models, areaSize, areaSize / 100);
	std::cout << "Iteration: " << i << "...\n";
    }

    std::cout << "Iterations count : " << stage << "\n";
    PrintData();
    Simulator::Destroy();

    return 0;
}
