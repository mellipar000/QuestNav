using System;
using NetworkTables;
using UnityEngine;

public class RayTracing : MonoBehaviour
{
    private Nt4Client _ntClient;

    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        // TODO: Change serverBaseAddress and potentially add serverPort, onOpen, and onNewTopicData
        _ntClient = new Nt4Client("Reefbreaker", "127.0.0.1");
        _ntClient.Connect();
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
