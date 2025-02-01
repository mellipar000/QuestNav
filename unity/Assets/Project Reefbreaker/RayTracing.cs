using System;
using NetworkTables;
using UnityEngine;

public class RayTracing : MonoBehaviour
{
    private Nt4Client _ntClient;
    // TODO: Input actual camera offset
    private const double cameraXOffset = 0.0, cameraYOffset = 0.0, cameraZOffset = 0.0, cameraPitch = 0.0;
    private double cameraXAngle = 0.0, cameraYAngle = 0.0, robotX = 0.0, robotY = 0.0, robotZ = 0.0;
    private String alliance = "Blue";
    private int _cameraAnglesSubscriptionId, _robotPositionSubscriptionId, _allianceSubscriptionId;

    // Field Constants
    private const double L1 = 0.4572, L2 = 0.70358, L3 = 1.106678, L4 = 1.652778;
    private readonly double[,] blueCoral = new double[,] {
        {3.822, 4.2}, {3.822, 3.859}, // A & B
        {4.0, 3.547}, {4.319, 3.372}, // C & D
        {4.66, 3.372}, {4.963, 3.557}, // E & F
        {5.138, 3.859}, {5.138, 4.2}, // G & H
        {4.963, 4.522}, {4.66, 4.678}, // I & J
        {4.319, 4.678}, {4.0, 4.522} // K & L
    };
    private readonly double[,] redCoral = new double[,] {
        {13.718, 3.859}, {13.718, 4.2}, // A & B
        {13.54, 4.522}, {13.221, 4.678}, // C & D
        {12.88, 4.678}, {12.577, 4.522}, // E & F
        {12.402, 4.2}, {12.402, 3.859}, // G & H
        {12.577, 3.557}, {12.88, 3.372}, // I & J
        {13.221, 3.372}, {13.54, 3.547} // K & L
    };
    private const double L23 = 0.905129, L34 = 1.308227;
    private readonly double[,] blueAlgae = new double[,] {
        {3.803, 4.025, L34}, // A & B
        {4.144, 3.42, L23}, // C & D
        {4.836, 3.42, L34}, // E & F
        {5.177, 4.025, L23}, // G & H
        {4.836, 4.61, L34}, // I & J
        {4.144, 4.61, L23} // K & L
    };
    private readonly double[,] redAlgae = new double[,] {
        {13.737, 4.025, L34}, // A & B
        {13.396, 4.61, L23}, // C & D
        {12.704, 4.61, L34}, // E & F
        {12.363, 4.025, L23}, // G & H
        {12.704, 3.42, L34}, // I & J
        {13.396, 3.42, L23} // K & L
    };

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        // Sets up network tables
        // TODO: Change serverBaseAddress and potentially add serverPort, onOpen, and onNewTopicData
        _ntClient = new Nt4Client("Reefbreaker", "127.0.0.1", onNewTopicData: OnNewTopicData);
        _ntClient.Connect();

        // Subscribes to relevant topics
        _cameraAnglesSubscriptionId = _ntClient.Subscribe("cameraAngles");
        _robotPositionSubscriptionId = _ntClient.Subscribe("robotPosition");
        _allianceSubscriptionId = _ntClient.Subscribe("alliance");
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    // Method run each time new values are posted to network tables
    public void OnNewTopicData(Nt4Topic topic, long timestamp, object value) {
        if (topic.Name == "cameraAngles" && value is double[] cameraAngles && cameraAngles.Length == 2) {
            cameraXAngle = cameraAngles[0];
            cameraYAngle = cameraAngles[1];
        } else if (topic.Name == "robotPosition" && value is double[] robotPose && robotPose.Length == 3) {
            robotX = robotPose[0];
            robotY = robotPose[1];
            robotZ = robotPose[2];
        } else if (topic.Name == "alliance" && value is String alliance) {
            this.alliance = alliance;
        }
    }
}
