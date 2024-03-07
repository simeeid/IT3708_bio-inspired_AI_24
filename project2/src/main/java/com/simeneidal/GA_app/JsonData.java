package com.simeneidal.GA_app;

import java.util.Map;

import com.google.gson.annotations.SerializedName;

public class JsonData {
    @SerializedName("instance_name")
    private String instanceName;
    @SerializedName("nbr_nurses")
    private int nbrNurses;
    @SerializedName("capacity_nurse")
    private int capacityNurse;
    private double benchmark;
    private Depot depot;
    private Map<String, Patient> patients;
    @SerializedName("travel_times")
    private double[][] travelTimes;

    // Getters and setters for the fields

    public String getInstanceName() {
        return instanceName;
    }

    public void setInstanceName(String instanceName) {
        this.instanceName = instanceName;
    }

    public int getNbrNurses() {
        return nbrNurses;
    }

    public void setNbrNurses(int nbrNurses) {
        this.nbrNurses = nbrNurses;
    }

    public int getCapacityNurse() {
        return capacityNurse;
    }

    public void setCapacityNurse(int capacityNurse) {
        this.capacityNurse = capacityNurse;
    }

    public double getBenchmark() {
        return benchmark;
    }

    public void setBenchmark(double benchmark) {
        this.benchmark = benchmark;
    }

    public Depot getDepot() {
        return depot;
    }

    public void setDepot(Depot depot) {
        this.depot = depot;
    }

    public Map<String, Patient> getPatients() {
        return patients;
    }

    public void setPatients(Map<String, Patient> patients) {
        this.patients = patients;
    }

    public double[][] getTravelTimes() {
        return travelTimes;
    }

    public void setTravelTimes(double[][] travelTimes) {
        this.travelTimes = travelTimes;
    }

    static class Depot {
        @SerializedName("return_time")
        private int returnTime;
        @SerializedName("x_coord")
        private int xCoord;
        @SerializedName("y_coord")
        private int yCoord;

        // Getters and setters for the fields
        public int getReturnTime() {
            return returnTime;
        }

        public void setReturnTime(int returnTime) {
            this.returnTime = returnTime;
        }

        public int getXCoord() {
            return xCoord;
        }

        public void setXCoord(int xCoord) {
            this.xCoord = xCoord;
        }

        public int getYCoord() {
            return yCoord;
        }

        public void setYCoord(int yCoord) {
            this.yCoord = yCoord;
        }
    }

    static class Patient {
        @SerializedName("x_coord")
        private int xCoord;
        @SerializedName("y_coord")
        private int yCoord;
        private int demand;
        @SerializedName("start_time")
        private int startTime;
        @SerializedName("end_time")
        private int endTime;
        @SerializedName("care_time")
        private int careTime;

        private int id;

        // Getters and setters for the fields

        public int getId() {
            return id;
        }

        public void setId(int id) {
            this.id = id;
        }

        public int getXCoord() {
            return xCoord;
        }

        public void setXCoord(int xCoord) {
            this.xCoord = xCoord;
        }

        public int getYCoord() {
            return yCoord;
        }

        public void setYCoord(int yCoord) {
            this.yCoord = yCoord;
        }

        public int getDemand() {
            return demand;
        }

        public void setDemand(int demand) {
            this.demand = demand;
        }

        public int getStartTime() {
            return startTime;
        }

        public void setStartTime(int startTime) {
            this.startTime = startTime;
        }

        public int getEndTime() {
            return endTime;
        }

        public void setEndTime(int endTime) {
            this.endTime = endTime;
        }

        public int getCareTime() {
            return careTime;
        }

        public void setCareTime(int careTime) {
            this.careTime = careTime;
        }
    }
}
