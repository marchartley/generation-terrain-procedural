#ifndef ENVOBJSERIALIZER_H
#define ENVOBJSERIALIZER_H

#include "EnvObject/EnvObject.h"
#include "EnvObject/EnvPoint.h"
#include "EnvObject/EnvCurve.h"
#include "EnvObject/EnvArea.h"
#include "EnvObject/EnvMaterial.h"

#include "Utils/json.h"

void to_json(nlohmann::json& json, const EnvObject& obj);
void to_json(nlohmann::json& json, const EnvPoint& obj);
void to_json(nlohmann::json& json, const EnvCurve& obj);
void to_json(nlohmann::json& json, const EnvArea& obj);
void to_json(nlohmann::json& json, const EnvMaterial& material);

void to_json(nlohmann::json& json, const DepositionRate& depos);
void to_json(nlohmann::json& json, const AbsorptionRate& absorb);


void from_json(const nlohmann::json& json, EnvObject& obj);
void from_json(const nlohmann::json& json, EnvPoint& obj);
void from_json(const nlohmann::json& json, EnvCurve& obj);
void from_json(const nlohmann::json& json, EnvArea& obj);
void from_json(const nlohmann::json& json, EnvMaterial& material);

void from_json(const nlohmann::json& json, DepositionRate& depos);
void from_json(const nlohmann::json& json, AbsorptionRate& absorb);

#endif // ENVOBJSERIALIZER_H
