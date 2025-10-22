# Script to clean spectrum file with proper decimal handling
input_filename = "Equipo1/Data/espectro1.txt"  # Change to your file path
output_filename = "Equipo1/Data/cleaned_spectrum.csv"

# Read the original file
with open(input_filename, 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find the data header line
data_start_line = -1
for i, line in enumerate(lines):
    if line.startswith("Pixel;Wavelength;Wavenumber;"):
        data_start_line = i
        header_line = line.strip()
        break

if data_start_line == -1:
    print("Error: Could not find data header line")
    exit()

# Extract column names from header
headers = header_line.split(';')
print("Available columns in file:")
for i, header in enumerate(headers):
    print(f"{i}: {header}")

# Ask user which columns to keep
wavelength_col = input("Enter the number for Wavelength column (default 1): ").strip()
intensity_col = input("Enter the number for Intensity column (default 6 for 'Raw data #1'): ").strip()

# Default to common positions if user just presses enter
wavelength_idx = int(wavelength_col) if wavelength_col else 1
intensity_idx = int(intensity_col) if intensity_col else 6

print(f"Using Wavelength column: {headers[wavelength_idx]}")
print(f"Using Intensity column: {headers[intensity_idx]}")

# Process data lines - fix the comma decimal issue
clean_data = []
for line in lines[data_start_line + 1:]:  # Start after header
    if line.strip() and not line.startswith(';') and not line.startswith('File'):
        parts = line.strip().split(';')
        if len(parts) > max(wavelength_idx, intensity_idx):
            wavelength = parts[wavelength_idx].strip().replace(',', '.')  # Replace comma with dot
            intensity = parts[intensity_idx].strip().replace(',', '.')    # Replace comma with dot
            
            # Only keep lines with valid data
            if wavelength and intensity and wavelength != '   ' and intensity != '   ':
                clean_data.append(f"{wavelength},{intensity}\n")

# Write clean data to new file
with open(output_filename, 'w', encoding='utf-8') as f:
    f.write("Wavelength,Intensity\n")  # New header
    f.writelines(clean_data)

print(f"Clean data saved to: {output_filename}")
print(f"Number of data points: {len(clean_data)}")

# Quick preview
print("\nFirst few lines of cleaned data:")
with open(output_filename, 'r', encoding='utf-8') as f:
    for i, line in enumerate(f):
        if i < 6:
            print(line.strip())
        if i == 5:
            print("...")
            break