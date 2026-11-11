CREATE SCHEMA IF NOT EXISTS derived;

DROP TABLE IF EXISTS derived.mock_outputevents CASCADE;
CREATE TABLE derived.mock_outputevents (
    subject_id integer,
    hadm_id integer,
    stay_id integer,
    charttime timestamp without time zone,
    storetime timestamp without time zone,
    itemid integer,
    value double precision,
    valueuom character varying(20)
);

DROP TABLE IF EXISTS derived.mock_inputevents CASCADE;
CREATE TABLE derived.mock_inputevents (
    subject_id integer,
    hadm_id integer,
    stay_id integer,
    starttime timestamp without time zone,
    endtime timestamp without time zone,
    storetime timestamp without time zone,
    itemid integer,
    amount double precision,
    amountuom character varying(20),
    rate double precision,
    rateuom character varying(20),
    orderid integer,
    linkorderid integer,
    patientweight double precision,
    ordercategorydescription character varying(250)
);

DROP TABLE IF EXISTS derived.mock_pharmacy CASCADE;
CREATE TABLE derived.mock_pharmacy (
    subject_id integer,
    hadm_id integer,
    pharmacy_id integer,
    poe_id character varying(25),
    starttime timestamp without time zone,
    stoptime timestamp without time zone,
    medication character varying(255),
    proc_type character varying(50),
    status character varying(50),
    enteral_route character varying(255),
    route character varying(255)
);
