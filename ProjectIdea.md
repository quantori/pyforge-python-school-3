# Project Idea

## Introduction

We were recommended by **Yaroslav Korkodinov** to add models:  Patients, Doctors and Drugs to the database. 

I was thinking about comming up with a project idea that would also include the molecules crud, that I already
spent much time on, and resemble a real world scenario.

So I came up with a hospital system, that is particularly focused on the personalization of the treatment,

Patients can be prescribed with the most suitable drugs and each drug can be composed of the most 
suitable molecules.

Maybe this does not sound very scientific and realistic, but I think it's at least nice for practicing. 
I had to learn about hierarchies in sqlalchemy and complex relationships. Additionally I used OAuth2, jwt, and role based access control.

## User Roles and Flows

### 1. Doctor

#### Registration Flow
1. **Fill Application Form:** Provide details such as email, password,  full name, experience, and other required information.
2. **Await Acceptance:** Wait until your account is reviewed and activated by hospital admins.
3. **Login:** Use the provided credentials to log into the system.

#### Seeing Patients
1. **Access Patient List:** View a list of all assigned patients.
2. **View Patient Details:** Access details for each patient (excluding  password).

#### Writing Prescriptions
1. **Select Patient:** Choose a patient for whom the prescription is to be written.
2. **Create Prescription:** Enter prescription details, including the drug (one drug per prescription) and other required fields.

### 2. Patient

#### Registration Flow
1. **Fill Application Form:** Provide details such as email, password, full name, and medical history.
2. **Await Acceptance:** Wait until your account is reviewed and activated by hospital admins.
3. **Login:** Use the provided credentials to log into the system.

#### Seeing Prescriptions
1. **Access Prescriptions:** View all prescriptions that have been prescribed to the patient.

#### Seeing Appointments
1. **View Appointments:** Access past and current appointments in the hospital.

#### Requesting Appointments
1. **Submit Request:** Provide details about the appointment request, including complaints, additional text, and preferred time.
2. **Await Confirmation:** Receive confirmation or scheduling details for the requested appointment.

### 3. Hospital Admin

#### Registration Flow
1. **Registration by Superadmin:** Superadmin is responsible for creating hospital admin accounts.

#### Seeing Registration Requests
1. **Access Requests:** View registration requests from patients and doctors.

#### Registering a Doctor/Patient
1. **Review Request:** Examine the registration request from a specific doctor or patient.
2. **Activate users:** Activate the accounts of registration requests

#### Managing Appointment Requests
1. **Access Requests:** View appointment requests from patients.
2. **Schedule Appointments:** Arrange appointments based on patient preferences and provided details.

### 4. Lab Admin

#### Managing Molecules and Drugs
1. **Populate System:** Enter and manage information related to molecules and drugs in the system.
2. **Substructure Search:** Perform searches for molecules based on substructures.

## Models

### DATABASE MODELS

![db_hospital - public.png](..%2F..%2FDesktop%2Fdb_hospital%20-%20public.png)

#### 1. Base 

- **created_at:** DateTime, default=now
- **updated_at:** DateTime, default=now, onupdate=now

#### 2. Molecule

- **id:** Integer, primary_key
- **name:** String, nullable=False
- **smiles:** String, nullable=False

#### 3. Drug

- **id:** Integer, primary_key
- **name:** String, nullable=False
- **description:** Text
- **molecules:** Many-to-Many relationship with Molecule

#### 4. User

- **id:** Integer, primary_key
- **email:** String, nullable=False
- **password:** String, nullable=False
- **full_name:** String, nullable=False
- **is_active:** Boolean, default=True
- **role:** Enum, values=['doctor', 'patient', 'hospital_admin', 'lab_admin']


#### 5. Patient(User)

- **gender:** Enum
- **date_of_birth:** Date
- **Medical_history:** Text

- **user_id:** Integer, ForeignKey('user.id')
- **doctor:** Many-to-One relationship with Doctor
- **prescriptions:** One-to-Many relationship with Prescription
- **appointments:** One-to-Many relationship with Appointment

#### 6. Doctor(User)

- **specialization:** String
- **experience:** Integer
- **description:** Text

- **user_id:** Integer, ForeignKey('user.id')
- **patients:** One-to-Many relationship with Patient
- **appointments:** One-to-Many relationship with Appointment


#### 7. appointment

- **doctor_id:** Integer, ForeignKey("doctor")
- **patient_id:** Integer, ForeignKey("patient")
- **date:** DateTime
- **complaints:** Text
- **additional_text:** Text
- **status:** Enum, values=['pending', 'confirmed', 'completed', 'cancelled']

#### 8. Prescription

- **doctor_id:** Integer, ForeignKey("doctor")
- **patient_id:** Integer, ForeignKey("patient")
- **drug_id:** Integer, ForeignKey("drug")
- **date:** DateTime
- **description:** Text


## API Endpoints

### `api/molecules/`

This endpoint is for lab admins. It allows them to manage molecules in the system.

Maybe some functionalities can be exposed to curious people, who want to know what 
substances are there in the medices they consume, or maybe to learn more about molecules.

#### **for lab_admin**:

    - POST: Create a new molecule

    - DELETE: /{molecule_id} Delete a molecule by id

    - PUT: /{molecule_id} Update a molecule by id

    - POST: /upload Upload a csv file containing molecules info


#### **for everyone**:
    
    - GET: /{molecule_id} Get a molecule by id

    - GET: /?page={page}&page_size={page_size}?name={name} Get all molecules, you can filter by name

    - GET: /search/smiles_containing_substructure?substructure={substructure_smiles} Get all molecules containing a given substructure

    - GET: /search/substructures_of?smiles={smiles} Get all substructures of a molecule

### `api/drugs/`
    
#### **for lab_admin**:

    - POST: Create a new drug

    - DELETE: /{drug_id} Delete a drug by id

    - PUT: /{drug_id} Update a drug by id

#### **for everyone**

    - GET: /?page={page}&page_size={page_size}?name={drug_name} Get all drugs, you can filter by name

    - GET: /{drug_id} Get a drug by id

### `api/patients/`
    
#### **for hospital_admin**:
    
    - GET: /{patient_id} Get a patient by id
    
    - POST: Register a new patient

    - DELETE: /{patient_id} Delete a patient by id

    - PUT: /{patient_id} Update a patient by id

#### **for patient**:
    
    - GET: /me get the patient's own profile

    - PATCH: /me update the patient's own profile info including password, but not email or role
    
    - GET: /me/prescriptions get all prescriptions of the patient

    - GET: /me/appointments?status={status}get all appointments of the patient you can filter by status, default is pending

    - POST: /me/appointments request a new appointment

### `api/doctors/`

#### **for hospital_admin**:
    
    - GET: /{doctor_id} Get a doctor by id
    
    - POST: Register a new doctor

    - DELETE: /{doctor_id} Delete a doctor by id

    - PUT: /{doctor_id} Update a doctor by id

#### **for doctor**:

    - GET: /me get the doctor's own profile

    - PATCH: /me update the doctor's own profile info including password, but not email or role
    
    - GET: /me/patients get all patients of the doctor

    - GET: /me/patients/{patient_id} get a patient of the doctor

    - POST: /me/patients/{patient_id}/prescriptions create a new prescription for a patient

    - GET: /me/appointments?status={status} get all appointments of the doctor, you can filter by status, default is pending


### `api/appointments/`

#### **for hospital_admin**:
    
    - GET: /{appointment_id} Get an appointment by id
    
    - GET: /?status={status} Get all appointments, you can filter by status

    - PUT: /{appointment_id} Update an appointment by id

    - DELETE: /{appointment_id} Delete an appointment by id
    

    



